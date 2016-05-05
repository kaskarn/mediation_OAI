library(magrittr); library(mice)
library(foreign); library(dplyr)
library(reshape2); library(ipw)
library(foreach)

whichtimes <- function (stub, data = df, plen = 3){
  grep(paste0(stub, "$"), colnames(data), value = TRUE) %>%
    substring(., 2, 3) %>% as.numeric(.)
}

#Exploration
CreateTableOne(data = df, vars = c("v00lfmaxf", "v00pase", "v00abcirc"), strata = "v00edcv")

#general function : pcs // adl
#knee function : womac disability // KOOS function
#physical exam : 20 m // isometric strength
#chronic disease : hypertension // diabetes

##### Read data from .csv downloaded at epi-ucsf.edu #####
tyears <- c(0,1,1.5,2,2.5,3:8)
times <- c(0:10)
df <- left_join(read.csv("Enrollees.txt", sep = "|"), read.csv("Outcomes99.txt", sep = "|") %>% rename(ID = id), by = "ID")
for(i in times) df <- left_join(df, read.csv(paste0("AllClinical", sprintf("%02d",i), ".txt"), sep = "|"), by = "ID")
colnames(df) <- tolower(colnames(df)) #fix SAS case-agnostic names
df %<>% subset(p02race %in% c(1,2)) %>% mutate(p02race = factor(p02race, labels = c("EA", "AA"))) %>%
  rename(v00bmi = p01bmi) #keep AA and EU

##### Define variables we'll keep for analyses #####
meta <- c("id", "version")
xsec <- c("p01sxkoa", "p02sex", "p02race", "v00edcv", "v00income")
long <- lapply(c("smkpkyr", "hsmss$", "bmi$", "20mpace$", "400mtim$", "^...age", "hspss$"), function (i) grep(i, colnames(df), value = TRUE))
outc <- c("v99elxioa", "v99erxioa", "v99rntcnt", "v99eddcf")
d2 <- df[,c(meta, xsec, unlist(long))]
for(i in colnames(d2)) if(is.integer(d2[[i]])) if(length(unique(d2[[i]])) < 10) d2[[i]] <- factor(d2[[i]])

##### infer missing age using time between T's #####
ageloc <- match("v00age", colnames(d2)) 
for(i in 1:length(times[-1])) d2[is.na(d2[ageloc+i]),ageloc+i] <- d2[is.na(d2[ageloc+i]),ageloc+i-1] + tyears[i+1] - tyears[i]

##### Meta stuff ##### 
N <- nrow(df)
impnum <- 5; bnum <- 50
outcomes <- c("20mpace", "400mtim", "hspss")
mediators <- c("obcat", "pycat", "depcat")
confounders <- c("age", "p02sex", "p02race", "agesq", "a_ra", "a2_ra", "a_se", "a2_se", "ra_se")
exposure <- "v00edcv"
len <- length(unique(d3[[exposure]]))

##### Multiple imputation using MICE #####
d2_imp <- mice(d2, m = impnum, pred=quickpred(d2,  include = c("v00age", "p02sex", "p02race"), exclude = meta)) %>%
  complete(., action = "long", include = FALSE)

##### convert to wide format #####
dl <- lapply(long, function (i) melt(d2_imp[match(c(".imp", meta, xsec, i), colnames(d2_imp))], 
                                     id.vars = c(meta, xsec, ".imp"), value.name = substring(i[1],4),
                                     variable.name = "time") %>% 
                                mutate(time = substr(i,2,3)[as.numeric(time)]))
d3 <- Reduce(function (...) merge(..., all = TRUE, by=c(meta, xsec, "time", ".imp")), dl)

##### Categorize mediators and create interaction terms##### 
d3 %<>% mutate(obcat = cut(bmi, breaks = c(min(bmi, na.rm = T),25,30,35, 999)),
              depcat = cut(hsmss, breaks = c(0, 40, 50, 60, 999)),
              pycat = .bincode(smkpkyr, breaks = c(0,10,20,30, 999)),
              agesq = age*age,
              ra_se = (p02race %in% "AA")*(p02sex == 2),
              a_ra = (p02race %in% "AA")*age,
              a2_ra = (p02race %in% "AA")*agesq,
              a_se = (p02sex == 2)*age,
              a2_se = (p02sex == 2)*agesq
              )
d3[d3$smkpkyr %in% 0,"pycat"] <- 0
d3$pycat <- factor(d3$pycat, labels = c("0", "0-10", "10-20", "20-30", "30+"))

# mediator formulas
med_models <- list(
  obcat = paste(exposure,"pycat", paste(confounders, collapse = "+"), sep = "+"),
  pycat = paste(exposure,"depcat", paste(confounders, collapse = "+"), sep = "+"),
  depcat = paste(exposure, paste(confounders, collapse = "+"), sep = "+")
)
med_formulas <- lapply(mediators, function (i) as.formula(paste(i, med_models[[i]], sep = "~")))

###### array storing unadjusted, controlled total effect, controlled direct effect overall ######
res_a <- array(NA, 
               dim=c(impnum, bnum, length(outcomes), 3, (len-1)), 
               dimnames = list(1:impnum, 1:bnum, outcomes, c("Unadjusted", "CTE", "CDE")))
###### array storing controlled direct effect for each mediator ######
res_b <- res_b_g <- array(NA, 
                    dim=c(impnum, bnum, length(outcomes), length(mediators), (len-1)), 
                    dimnames = list(1:impnum,1:bnum, outcomes, mediators))

###### Analyses education #####

exposure <- "v00edcv"
pb <- txtProgressBar(style=3); fin <- 0
d4 <- d3[d3$.imp %in% 1:3,]
for(i in levels(d3$.imp)){
  #timp <- d3[d3$.imp == i,]
  for(j in 1:bnum){
    tsamp <- d3[d3$.imp == i,][d3[d3$.imp == i,]$id %in% sample(df$id, N, replace = TRUE),]
    wa_num <-  multinom(get(exposure) ~ 1, data = tsamp, trace = FALSE) %$% 
      predict(object = ., newdata = tsamp, "probs")
    wa_den <- multinom(get(exposure) ~ p02sex + p02race + age + age*age + p02sex*p02race + 
                         age*p02sex + age*p02race, data = tsamp, trace = FALSE) %$% 
              predict(object = ., newdata = tsamp, "probs")
    tsamp$wa <- (wa_num / wa_den)[cbind(1:nrow(tsamp),as.numeric(tsamp[[exposure]]))]
    
    wm_num.models <- lapply(mediators, function (i) multinom(data = tsamp, get(i) ~ get(exposure), trace = FALSE))
    wm_den.models <- lapply(1:length(mediators), function (i) 
      multinom(data = tsamp, med_formulas[[i]], trace = FALSE))
    for(m in 1:length(mediators)){
      wm_tmp <- predict(wm_num.models[[m]], tsamp, "probs")/predict(wm_den.models[[m]], tsamp, "probs")
      tsamp[[paste0("wm_", mediators[m])]] <- wm_tmp[cbind(1:nrow(tsamp),as.numeric(tsamp[[mediators[m]]]))]
    }
    for(o in 1:length(outcomes)){
      for(k in len){
        
      }
      res_a[i,j,o,1,] <- lm(data = tsamp, get(outcomes[o]) ~ get(exposure))$coefficients[2:(len)]
      res_a[i,j,o,2,] <- lm(data = tsamp, get(outcomes[o]) ~ get(exposure), weights = wa)$coefficients[2:len]
      for(m in 1:length(mediators)){
        res_b[i,j,o,m,] <- lm(data = tsamp, get(outcomes[o]) ~ get(exposure) + get(mediators[m]), 
                              weights = wa * get(paste0("wm_", mediators[m])))$coefficients[2:(len)]
      }
    }
    fin <- fin + 1
    setTxtProgressBar(pb, fin/length(levels(d3$.imp))/bnum)
  }
}

##### present results #####
for(i in 1:length(mediators)){
  print(mediators[i])
  for(o in 1:length(outcomes)){
    print(outcomes[o])
    print(apply((res_a[,,o,2,] - res_b[,,o,i,])/res_a[,,o,i,], 3, quantile, probs = c(0.025, 0.5, 0.975))) #
  }
}


#functionalizing
solve_graphA <- function(df, x, y, m, c, z, imp = ".imp", id = "id", N, nboot = 100){
  impnum <- length(unique(df[[imp]]))
  res <- array(NA, dim = c(impnum, nboot, 3, length(levels(df[[x]]))))
  for(i in levels(df[[imp]])){
    for(j in nboot){
      tsamp <- df[df[[imp]] == i,][df[df[[imp]] == i,][[id]] %in% sample(unique(df[[id]]), N, replace = TRUE)]
      wa_num <-  multinom(get(exposure) ~ 1, data = tsamp, trace = FALSE) %$% 
        predict(object = ., newdata = tsamp, "probs")
      wa_den <- multinom(as.formula(paste0(exposure, "~", paste0(c, collapse = "+"))), data = tsamp, trace = FALSE) %$% 
        predict(object = ., newdata = tsamp, "probs")
      tsamp$wa <- (wa_num / wa_den)[cbind(1:nrow(tsamp),as.numeric(tsamp[[exposure]]))]
    }
  }
}

# 
# ### g-formula stuff ###
# tm_med <- tsamp[[mediators[m]]]; tm_x <- tsamp[[exposure]]
# tsamp[[exposure]] <- unique(tsamp$v00edcv)[1]
# tsamp[[mediators[m]]] <- predict(wm_den.models[[m]], tsamp, "class")
# tsamp[[exposure]] <- tm_x
# tm_out <- lm(data = tsamp, get(outcomes[o]) ~ get(exposure) + get(mediators[m]) + p02sex + p02race + age + age*age) %$%
#   predict(., tsamp)
# res_b_g[i,j,o,m,] <- mean(tm_out - tsamp[[outcomes[o]]])
