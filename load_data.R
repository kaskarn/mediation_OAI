library(magrittr); library(mice)
library(foreign); library(dplyr)
library(reshape2); library(ipw)
library(foreach); library(nnet)

## Summarizes data availability
whichtimes <- function (stub, data = df, plen = 3){
  grep(paste0(stub, "$"), colnames(data), value = TRUE) %>%
    substring(., 2, 3) %>% as.numeric(.)
}

# Exploration
#CreateTableOne(data = df, vars = c("v00lfmaxf", "v00pase", "v00abcirc"), strata = "v00edcv")

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
