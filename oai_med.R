library(magrittr); library(mice)
library(foreign); library(dplyr)
library(reshape2); library(ipw)
library(foreach); library(nnet)

X <- "v00edcv" #categorical education
xlen <- length(levels(dat[[X]]))-1
M <- "obcat" # obesity categories. 
Y <- "hspss"
L <- "smoker" #smoking status
C <- c("age", "agesq", "p02race", "p02sex", "ra_se")
boot <- 10
res <- array(NA, dim = c(boot, xlen, 3))

#### RANDOM INTERVENTIONAL ANALOGUE ######

pbar <- txtProgressBar(style = 3)
for(i in 1:boot){
  dat <- sample_frac(filter(xs, p01sxkoa != 0), 1, replace = TRUE)
  while(min(table(dat[[X]])) < 50) dat <- sample_frac(xs, 1, replace = TRUE)
  
  rint_items <- rint_med.mkdata(dat, X = X, C = C, M = M, L = L)
  dec <- rint_med.decompose(rint_items$an_dat, Y = Y, X = X)
  res[i,,1] <- dec$nder
  res[i,,2] <- dec$nier
  res[i,,3] <- dec$ter
  setTxtProgressBar(pbar, i/boot)
}

res95 <- lapply(1:3, function (i) apply(res[,,i], 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm =TRUE))




##### WEIGHTING-BASED METHOD FOR CDE #####

test_0 <- med_iptw.mkdat(filter(xs, p01sxkoa == 0), A = X, M = M, Y = Y, C = C, L = L)
test_1 <- med_iptw.mkdat(filter(xs, p01sxkoa != 0), A = X, M = M, Y = Y, C = C, L = L)
test_a <- med_iptw.mkdat(xs, A = X, M = M, Y = Y, C = C, L = L)

#res <- med_iptw.decomp(test, A = X, M = M, Y = Y)
res_0_noint <- rbind(
  med_iptw.decomp(test_0, A = X, M = M, Y = Y, noint = TRUE),
  lm(data = test_0, hspss ~ v00edcv, weights = ipw_conf)$coefficients[-1]
)
res_1_noint <- rbind(
  med_iptw.decomp(test_1, A = X, M = M, Y = Y, noint = TRUE),
  lm(data = test_0, hspss ~ v00edcv, weights = ipw_conf)$coefficients[-1]
)
res_a_noint <- rbind(
  med_iptw.decomp(test_a, A = X, M = M, Y = Y, noint = TRUE),
  lm(data = test_0, hspss ~ v00edcv, weights = ipw_conf)$coefficients[-1]
)



##### Regression-based #####
# X <- "bmi"
# M <- "hsmss"
# C <- confounders
# 
# test <- twoway.cont.cont(xs, X = "exposure", M = "bmi", Y = "hspss", C = confounders, mlvl = 25:30, noint = FALSE)
# test <- twoway.cont.cont(xs, X = "exposure", M = "bmi", Y = "hspss", C = confounders, mlvl = 25:30, noint = FALSE)