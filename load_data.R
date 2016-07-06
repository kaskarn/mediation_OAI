library(magrittr)
library(foreign)
library(dplyr)

## Download files


##### Read ALL data from .csv files downloaded at epi-ucsf.edu (allclinical, outcomes, enrollees). Version variables gets@
#messed up but that does not matter
tyears <- c(0,1,1.5,2,2.5,3:8)
times <- c(0:10)
df <- left_join(read.csv("data/Enrollees.txt", sep = "|"), 
                read.csv("data/Outcomes99.txt", sep = "|") %>% rename(ID = id), by = "ID")

for(i in times) df <- left_join(df, read.csv(paste0("data/AllClinical", sprintf("%02d",i), ".txt"), sep = "|"), by = "ID")
colnames(df) <- tolower(colnames(df)) #fix SAS case-agnostic names


## Fix 20m and 400m variable suffix to start with letter. issue is with SAS/R naming conventions
## New convention is v##w20mspace and v##w400mtim
names(df)[grep("20mpace", names(df))] <- unlist(lapply(grep("20mpace", names(df), value = TRUE), 
                                                       function (i) paste0(substring(i,1,3), "w",substring(i,4))))

names(df)[grep("400mtim", names(df))] <- unlist(lapply(grep("400mtim", names(df), value = TRUE), 
                                                       function (i) paste0(substring(i,1,3), "w",substring(i,4))))

## Define variables we'll keep in the reduced dataset, grabbing longitudinal variables from their suffix with grep() 
meta <- c("id", "version")
xsec <- c("p01sxkoa", "p02sex", "p02race", "v00edcv", "v00income", "p01xrkoa", "p02elgrisk", "p02wtga", "p02cncr3")
food <- c("ffq66$", "dtna$", "dtkcal", "srvfrt", "dtchol")
long <- lapply(c("hsmss$", "bmi$", "w20mpace$", "w400mtim$", "^...age", "hspss$", "^...smoker$", "smkpkyr",
                 "bpsys$", "bpdias$", food, "abcirc$", "smkamt$", "drnkamt$"), function (i) grep(i, colnames(df), value = TRUE))
outc <- c("v99elxioa", "v99erxioa", "v99rntcnt", "v99eddcf")

## Build longitudinal dataset
for(i in unlist(c(meta, xsec, long))) if(!(i %in% names(df))) stop(paste(i, "not a variable in dataset"))
d2 <- df[c(meta, xsec, unlist(long))]

### Basic data fixes (not new variables)
## Data fixes and reformats
# All discrete ints to factor
for(i in colnames(d2)) if(is.integer(d2[[i]])) if(length(unique(d2[[i]])) < 10) d2[[i]] <- factor(d2[[i]])

# Smoking categories:
for(i in grep("smoker$", colnames(d2))) levels(d2[[i]]) <- c("Never", "Current", "Former", "Never") #Second "Never" is occasional, not regular ever smokers)
for(i in grep("drnkamt$", colnames(d2))) levels(d2[[i]]) <- c("None", "1-3", "4-7", "8-14", "15-21", "22-27", "28+")

# Infer missing age using time between T's:
ageloc <- match("v00age", colnames(d2))
for(i in 1:length(times[-1])) d2[is.na(d2[ageloc+i]),ageloc+i] <- d2[is.na(d2[ageloc+i]),ageloc+i-1] + tyears[i+1] - tyears[i]

# Categorize education and change referent
d2$v00edcv %<>% factor(labels = c("<HS", "HS", "S.C.", "Coll", "S.G.", "Grad"))
d2$v00edcv %<>% relevel(ref = last(levels(.)))

## Save longitudinal dataset
write.csv(d2, "oai_long_july.csv")
###############################

#### Create cross-sectiona data 
## Only keep first timepoint of measurements and drop prefix
long_xs <- sapply(long, "[[", 1)
xs <- d2[c(meta, xsec, long_xs)]
names(xs)[match(long_xs, names(xs))] %<>% substring(., 4)

## Add new variables; these will need to be passively imputed in longitudinal dataset
xs %<>% mutate( obcat = cut(bmi, breaks = c(min(bmi, na.rm = T), 20,25,30,35, max(bmi, na.rm = T))) %>% factor,
                depcat = cut(hsmss, breaks = c(0, 40, 50, 60, 999)) %>% factor,
                pycat = .bincode(smkpkyr, breaks = c(0,10,20,30, 999)) %>% factor,
                agesq = age*age,
                ra_se = and((p02race %in% "AA"),(p02sex == 2)),
                a_ra = (p02race %in% "AA")*age,
                a2_ra = as.numeric(p02race %in% "AA")*agesq,
                a_se = as.numeric(p02sex == 2)*age,
                a2_se = as.numeric(p02sex == 2)*agesq
              )
for(i in c("obcat", "depcat", "pycat")) xs[[i]] <- factor(xs[[i]])

## Set referent obesity category to 'healthy' level (20-25)
xs$obcat <- relevel(xs$obcat, ref = levels(xs$obcat)[2])

## Crate additional smoking variables dummy variables (current, ever)
xs$smknow <- xs$smkever <- xs$smoker
levels(xs$smknow) <- c("Never", "Current", "Never", "Never")
levels(xs$smkever) <- c("Never", "Ever", "Ever", "Never")

## Save cross-sectional dataset
write.csv(xs, "oai_xs_july.csv")
###############################
