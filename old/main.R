

#load data and format
source("load_data.R")

##### Meta 
N <- nrow(df)
impnum <- 1; bnum <- 20
outcomes <- c("hspss", "w20mpace") #"20mpace", "400mtim", 
mediators <- c("bmi", "hsmss")
confounders <- c("age", "p02sex", "p02race", "agesq", "a_ra", "a2_ra", "a_se", "a2_se", "ra_se")
exposure <- "v00edcv"
xlen <- levels(d2[[exposure]]) %>% length - 1

set.seed(0808)


########## TESTING AREA ############### 

## randomized interventional analogue testing

X <- "v00edcv" #categorical education
M <- "obcat" # obesity categories. 
Y <- "hspss"
L <- "smoker" #smoking status

#test!
test <- rint_med(xs, X = "v00edcv", M = "obcat", Y = "hspss", C = confounders, L = "smoker")

#######################################

# 
# ### tests of regresson methdos
# 

# ########################################   