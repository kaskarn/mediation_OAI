library(magrittr); library(mice)
library(foreign); library(dplyr)
library(reshape2); library(ipw)
library(foreach); library(nnet)

## Summarizes data availability (for exploration)
whichtimes <- function (stub, data = df, plen = 3){
  grep(paste0(stub, "$"), colnames(data), value = TRUE) %>%
    substring(., 2, 3) %>% as.numeric(.)
}
fmatch <- function(v, x, l = -1, mat = FALSE){
  #replaces syntax of type [cbind(1:nrow(dat),match(dat[[M]], levels(dat[[M]])))]
    if(!is.factor(x)) stop("Must be factor variable")
    if(l == -1) l <- length(x)
    
    if(length(levels(x)) == 2) v <- cbind(1-v, v)
  
    if(mat){return(cbind(1:l,match(x, levels(x))))
    }else return(v[cbind(1:l,match(x, levels(x)))])
}
## Twoway decomposition indirect and direct effects  (VDW)
 #for continuous mediator and outcome

twoway.cont.cont <- function(df, X, M, Y, C = "", xlen = 0, mlen = 0, noint = FALSE, mlvl = NULL, delta = FALSE){
  #sets position and lengths of coefficient vectors
  if(xlen == 0) if(is.factor(df[[X]])) xlen <- levels(df[[X]]) %>% length - 1 else xlen <- 1
  if(mlen == 0) if(is.factor(df[[M]])) mlen <- levels(df[[M]]) %>% length - 1 else mlen <- 1
  c_mu <- lapply(C, function (i) if(is.factor(df[[i]])) table(df[[i]])[-1]/nrow(df) else mean(df[[i]])) %>% unlist
  
  #throws error if needed
  if((xlen > 1 | mlen > 1) & delta == TRUE) stop("Standard error must be estimated by bootstrapping for polytomous treatment or mediator") 
  
  #adds parentheses to allow flexible model specification
  for(i in c("X","M","Y","C")) assign(i, paste0("(",get(i),")"))
  
  #model of mediator
  m1 <- lm(data = df, as.formula(paste(M, paste(X, paste(C, collapse = "+"), sep = "+"),  sep = "~")))
  
  #model of outcome
  if(noint) { m2 <- lm(data = df, as.formula(paste(Y, paste(X, M, paste(C, collapse = "+"), sep = "+"),  sep = "~")))
  }else m2 <-  lm(data = df, as.formula(paste(Y, paste(X, M, paste(C, collapse = "+"), paste(X,M, sep = "*"), sep = "+"),  sep = "~")))
  
  #breaks up results into different variable coefficients
  p1 <- cumsum(c(1,xlen,length(c_mu)))
  p2 <- cumsum(c(1,xlen,mlen,length(c_mu),mlen*xlen))
  
  #creates shorthands for coefficients to be used
  b0 <- m1$coefficients[1]
  for(i in 1:4) assign(paste0("t",i), c(m2$coefficients, rep(0,xlen*mlen*noint))[(p2[i]+1):p2[i+1]] %>% unname)
  for(i in 1:2) assign(paste0("b",i), m1$coefficients[(p1[i]+1):p1[i+1]] %>% unname)
  
  #output results, optionally with controlled direct effects
  out <- list(
    nde = (t1 + t4*b0 + t4*b1*ref + t4*b2%*%c_mu)*(treat - ref),
    nie = (t2*b1 + t4*b1*treat)*(treat - ref),
    #throws mediated moderation in there, although not two-way decomposition
    mint = (t4*b1)*(treat - ref)
  )
  if (!is.null(mlvl)) out <- list(out, cde = do.call(rbind, lapply(mlvl, function (i) ((t1 + t4*i)*(treat - ref)) )) )
  
  #Delta method
  if(delta = TRUE){
    sigma <- cbind(
      rbind(vcov(m1), matrix(0, length(m2$coefficients), length(m1$coefficients))),
      rbind(vcov(m2), matrix(0, length(m1$coefficients), length(m2$coefficients)))
    )
    #same for Gamma
    gamma <- list(  
      nde = c(t4, t4*ref, t4*c_mu, 0, 1, 0, b0 + b1*ref + b2%*%c_mu, rep(0, length(c_mu))),
      nie = c(0, t2 + t4*treat, rep(0, length(c_mu)), 0, 0, b1, b1*treat, rep(0, length(c_mu)))
    )
    if (!is.null(mlvl)) gamma <- list(gamma, cde = do.call(rbind, lapply(mlvl, function (i) c(rep(0, length(c_mu)+3), 1, 0, i, rep(0, length(c_mu))))))
    se <- lapply(c("nde", "nie"), function (i) sqrt(gamma[[i]]%*%sigma%*%gamma[[i]]) * abs(treat - ref))
  }
  
  #end
  return(list(out = out, se = se))
}

rint_w <- function(dat, X, M, Y, C = "", L, xlen = 0, mlen = 0, ref = 0, treat = 1, stable = TRUE){
  ###Create referent indicator A = a*
  dat$aref <- as.integer(dat[[X]] == levels(dat[[X]])[1])
  ref <- levels(dat[[X]])[1]

  ###replicate input with parentheses
  for(i in c("X","M","Y","C", "L")) assign(paste0(i, "2"), paste0("(",get(i),")"))
  
  ### Compute weights
  ## denominator
  # p(a | c)
  a_form <- as.formula(paste0(X, "~", paste0(C2, collapse = "+")))
  a_mod <- multinom(data = dat, a_form, trace = FALSE)
  a_pro <- predict(object = a_mod, newdata = dat, "probs") %>% fmatch(dat[[X]])
  
  # p(m |l, a, c)
  m_form <- as.formula(paste0(M, "~", paste(X, paste0(L2, collapse = "+"), paste0(C2, collapse = "+"), sep = "+")))
  m_mod <- multinom(data = dat, m_form1, trace = FALSE)
  m_pro <- predict(object = m_mod1, newdata = dat, "probs") %>% fmatch(dat[[M]])

  den <- a_pro * m_pro
  
  ## numerator
  # p(l | a*, c)
  l_form <- as.formula(paste0(L2, "~", X, "+", paste0(C2, collapse = "+")))
  l_mod <- multinom(data = dat, l_form, trace = FALSE)
  l_pro <- predict(object = l_mod, newdata = mutate_(dat, .dots = setNames(list(~ref), X)), "probs")
  
  # math to get numerator
  num <- lapply(levels(dat[[L]]), 
                function (i) {
                    predict(
                      object = m_mod, 
                      newdata = mutate_(dat, .dots = setNames(list(~i), L)) %>% 
                                mutate_(., .dots = setNames(list(~ref), X)),
                      "probs"
                    ) %>% fmatch(dat[[M]]) * l_pro[,match(i, levels(dat[[L]]))]
                }) %>% do.call(cbind, .) %>% apply(., 1, sum)
  
  #get weight
  if(stable){ dat$w <- num/den * ((table(dat[[X]])/nrow(dat))[match(dat[[X]], levels(dat[[X]]))])
  }else dat$w <- num/den  

  ###duplicate dataset for each category of exposure
  dat <- do.call(rbind, lapply(1:length(levels(d2[[X]])), function(i) mutate(dat, astar = i - 1))) %>% 
         mutate(astar = factor(astar, labels = levels(dat[[X]])))
  
#   for(i in levels(dat[[X]])) dat[[paste0("astar_",match(i, levels(dat[[X]])))]] <- i == dat$astar
#   asl <- grep("^astar_", names(dat), value = TRUE)
  
#   N <- nrow(dat)
#   dat2 <- dat
#   dat2$astar <- dat[[X]]
#   dat$astar = NA
#   dat <- rbind(dat2, dat)
#   dat[is.na(dat$astar),"astar"] <- levels(dat[[X]])[1]
  
  
  ### direct effect
  nde_mod <- lm(data = dat[dat$astar == ref,], get(Y) ~ get(X), weights = w)
  nie_mod <-lm(data = dat[dat[[X]] != ref,], get(Y) ~ astar, weights = w)
}

X <- exposure
X <- "ed2"
M <- "obcat"
Y <- "hspss"
L <- "smoker"
test <- twoway.cont.cont(xs, X = "exposure", M = "bmi", Y = "hspss", C = confounders, mlvl = 25:30, noint = FALSE)

test <- twoway.cont.cont(xs, X = "exposure", M = "bmi", Y = "hspss", C = confounders, mlvl = 25:30, noint = FALSE)

X <- "bmi"
M <- "hsmss"
C <- confounders

## threeway decomposition
# output is vector of length equal to number of categories of outcome, minus 1 

##### Read data from .csv downloaded at epi-ucsf.edu #####
tyears <- c(0,1,1.5,2,2.5,3:8)
times <- c(0:10)
df <- left_join(read.csv("Enrollees.txt", sep = "|"), read.csv("Outcomes99.txt", sep = "|") %>% rename(ID = id), by = "ID")
for(i in times) df <- left_join(df, read.csv(paste0("AllClinical", sprintf("%02d",i), ".txt"), sep = "|"), by = "ID")
colnames(df) <- tolower(colnames(df)) #fix SAS case-agnostic names
df %<>% subset(p02race %in% c(1,2)) %>% mutate(p02race = factor(p02race, labels = c("EA", "AA"))) %>%
  rename(v00bmi = p01bmi) #keep AA and EU %>% 
#fix 20m and 400m variables suffix to start with letter
names(df)[grep("20mpace", names(df))] <- unlist(lapply(grep("20mpace", names(df), value = TRUE), function (i) paste0(substring(i,1,3), "w",substring(i,4))))
names(df)[grep("400mtim", names(df))] <- unlist(lapply(grep("400mtim", names(df), value = TRUE), function (i) paste0(substring(i,1,3), "w",substring(i,4))))

##### Define variables we'll keep for analyses #####
meta <- c("id", "version")
xsec <- c("p01sxkoa", "p02sex", "p02race", "v00edcv", "v00income")
long <- lapply(c("hsmss$", "bmi$", "w20mpace$", "w400mtim$", "^...age", "hspss$", "^...smoker$", "smkpkyr"), function (i) grep(i, colnames(df), value = TRUE))
outc <- c("v99elxioa", "v99erxioa", "v99rntcnt", "v99eddcf")
d2 <- df[c(meta, xsec, unlist(long))][df$p01sxkoa != 0,]

##### Data fixes and reformats
# all discrete ints to factor
for(i in colnames(d2)) if(is.integer(d2[[i]])) if(length(unique(d2[[i]])) < 10) d2[[i]] <- factor(d2[[i]]) 
#smoking categories:
for(i in grep("smoker$", colnames(d2))) levels(d2[[i]]) <- c("Never", "Current", "Former", "Never") #Second "Never" is occasional, not regular ever smokers)
#infer missing age using time between T's: 
ageloc <- match("v00age", colnames(d2)) 
for(i in 1:length(times[-1])) d2[is.na(d2[ageloc+i]),ageloc+i] <- d2[is.na(d2[ageloc+i]),ageloc+i-1] + tyears[i+1] - tyears[i]
#categorize education and change referent
d2$v00edcv %<>% factor(., labels = c("<HS", "HS", "S.C.", "Coll", "S.G.", "Grad"))
d2$v00edcv  %<>% relevel(ref = last(levels(.)))

#create baseline data
long_xs <- sapply(long, "[[", 1)
xs <- d2[c(meta, xsec, long_xs)] 
names(xs)[match(long_xs, names(xs))] %<>% substring(., 4)
xs %<>% mutate(obcat = cut(bmi, breaks = c(min(bmi, na.rm = T),25,30,35, 999)),
                depcat = cut(hsmss, breaks = c(0, 40, 50, 60, 999)),
                pycat = .bincode(smkpkyr, breaks = c(0,10,20,30, 999)),
                agesq = age*age,
                ra_se = (p02race %in% "AA")*(p02sex == 2),
                a_ra = (p02race %in% "AA")*age,
                a2_ra = (p02race %in% "AA")*agesq,
                a_se = (p02sex == 2)*age,
                a2_se = (p02sex == 2)*agesq
              )
xs$smknow <- xs$smkever <- xs$smoker 
levels(xs$smknow) <- c("Never", "Current", "Never", "Never")
levels(xs$smkever) <- c("Never", "Ever", "Ever", "Never")

##### Meta for analyses
N <- nrow(df)
impnum <- 1; bnum <- 20
outcomes <- c("hspss", "w20mpace") #"20mpace", "400mtim", 
mediators <- c("bmi", "hsmss")
confounders <- c("age", "p02sex", "p02race", "agesq", "a_ra", "a2_ra", "a_se", "a2_se", "ra_se")
exposure <- "v00edcv"
xlen <- levels(d2[[exposure]]) %>% length - 1

set.seed(0808)

##### Cross sectional analyses #####

res_int_xs <- res_noint_xs <- list(outcomes)
te <- nie <- nde <- matrix(nrow = bnum, ncol = xlen)
colnames(te) <- colnames(nie) <- colnames(nde) <- levels(xs[[exposure]])[-1]
for(o in outcomes){
  res_int_xs[[o]] <- list(mediators)
  for(m in mediators){
    for(i in 1:bnum){
      tsamp <- sample_frac(xs, 1, replace = TRUE)
      te[i,] <-  lm(data = tsamp, as.formula(paste0(o, "~", paste(exposure, paste(confounders, collapse = "+"), sep = "+"))))$coefficients[2:(xlen+1)]
      nie[i,] <- indirect.reg.cont(df = tsamp, X = exposure, M = m, Y = o, C = confounders, xlen = xlen, mlen = 1)
      nde[i,] <- te[i,] - nie[i,]
    }
    res <- lapply(c("te", "nie", "nde"), function (i) list(crude = get(i), result = apply(get(i), 2, quantile, probs = c(0.025, 0.50, 0.975))))
    res_int_xs[[o]][[m]] <- list(te = res[[1]], nie = res[[2]], nde = res[[3]])
  } 
} 
#####
##### Longitudinal analyses ##### 
meta <- c(meta, "tok")
pbar <- txtProgressBar(style=3)
res_int <- array(NA, dim = c(bnum, length(outcomes), 3, xlen))
res_noint <- array(NA, dim = c(bnum, length(outcomes), 3, xlen))

for(i in 1:bnum){
  ## Resample
  tsamp <- sample_frac(d2, 1, replace = TRUE) %>% group_by(id) %>% mutate(tok = row_number())
  
  ## Impute
  t_imp <- mice(tsamp, m = 1, maxit = 5, 
                pred=quickpred(tsamp,  include = c("v00age", "p02sex", "p02race"), exclude = meta)) %>%
                complete(., action = "long", include = FALSE)
  ## Convert to long
  dl <- lapply(long, function (j) melt(t_imp[match(c(meta, xsec, j), colnames(t_imp))], 
                                       id.vars = c(meta, xsec), value.name = substring(j[1],4),
                                       variable.name = "time") %>% 
                                  mutate(time = substr(j,2,3)[as.numeric(time)],
                                         time = as.numeric(time)))

  tsamp <- Reduce(function (...) merge(..., all = TRUE), dl)
  
  ## Categorize mediators and create interaction terms
  tsamp %<>% mutate(obcat = cut(bmi, breaks = c(min(bmi, na.rm = T),25,30,35, 999)),
                    depcat = cut(hsmss, breaks = c(0, 40, 50, 60, 999)),
                    pycat = .bincode(smkpkyr, breaks = c(0,10,20,30, 999)),
                    agesq = age*age,
                    ra_se = (p02race %in% "AA")*(p02sex == 2),
                    a_ra = (p02race %in% "AA")*age,
                    a2_ra = (p02race %in% "AA")*agesq,
                    a_se = (p02sex == 2)*age,
                    a2_se = (p02sex == 2)*agesq
                   )
  tsamp[tsamp$smkpkyr %in% 0,"pycat"] <- 0
  tsamp$exp_stor <- tsamp[[exposure]] 
  tsamp$pycat <- factor(tsamp$pycat, labels = c("0", "0-10", "10-20", "20-30", "30+"))

  for(o in 1:length(outcomes)){
    res_int[i, o, 1,] <- lm(data = tsamp, as.formula(paste0(outcomes[o], "~", paste(exposure, paste(confounders, collapse = "+"), sep = "+"))))$coefficients[2:(xlen+1)]
    res_int[i, o, 3,] <- indirect.reg.cont(df = tsamp, X = exposure, M = "bmi", Y = outcomes[o], C = confounders, xlen = xlen, mlen = 1)
    res_noint[i, o, 1,] <- res_reg_int[i, o, 1,]
    res_noint[i, o, 3,] <- indirect.reg.cont(df = tsamp, X = exposure, M = "bmi", Y = outcomes[o], C = confounders, xlen = xlen, mlen = 1, noint = TRUE)
  }
  setTxtProgressBar(pbar, i/bnum)
}
