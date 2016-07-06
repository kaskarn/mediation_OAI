library(magrittr)
library(nnet)
library(dplyr)
boot <- 1000
#levels(xs$obcat) <- levels(xs$obcat)[c(1,1,3,4,5)]

res_multi <- list(
  roa = med_multi(filter(xs, p01xrkoa != 0), Y = hspss, A = v00edcv, M = obcat + smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  soa = med_multi(filter(xs, p01sxkoa != 0), Y = hspss, A = v00edcv, M = obcat + smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  noa = med_multi(filter(xs, p01xrkoa == 0), Y = hspss, A = v00edcv, M = obcat + smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  all = med_multi(xs, Y = hspss,  A = v00edcv, M = obcat + smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot)
)

res_bmi <- list(
  roa = med_iptw(filter(xs, p01xrkoa != 0), Y = hspss, A = v00edcv, M = obcat, L = smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  soa = med_iptw(filter(xs, p01sxkoa != 0), Y = hspss, A = v00edcv, M = obcat, L = smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  noa = med_iptw(filter(xs, p01xrkoa == 0), Y = hspss, A = v00edcv, M = obcat, L = smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  all = med_iptw(xs, Y = hspss, A = v00edcv, M = obcat, L = smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot)
)

res_bmi2 <- list(
  roa = med_smean(filter(xs, p01xrkoa != 0), Y = hspss, A = v00edcv, M = bmi, L = smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  soa = med_smean(filter(xs, p01sxkoa != 0), Y = hspss, A = v00edcv, M = bmi, L = smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  noa = med_smean(filter(xs, p01xrkoa == 0), Y = hspss, A = v00edcv, M = bmi, L = smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  all = med_smean(xs, Y = hspss, A = v00edcv, M = bmi, L = smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot)
)
  
res_bmi2 <- list(
  roa = med_smean(filter(xs, p01xrkoa != 0), Y = hspss, A = v00edcv, M = abcirc, L = smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  soa = med_smean(filter(xs, p01sxkoa != 0), Y = hspss, A = v00edcv, M = abcirc, L = smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  noa = med_smean(filter(xs, p01xrkoa == 0), Y = hspss, A = v00edcv, M = abcirc, L = smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  all = med_smean(xs, Y = hspss, A = v00edcv, M = abcirc, L = smoker + dtna + drnkamt, C = age + agesq + p02sex + p02race, boot = boot)
)

res_smk <- list(
  roa = med_iptw(filter(xs, p01xrkoa != 0), Y = hspss, A = v00edcv, M = smoker, C = age + agesq + p02sex + p02race, boot = boot),
  soa = med_iptw(filter(xs, p01sxkoa != 0), Y = hspss, A = v00edcv, M = smoker, C = age + agesq + p02sex + p02race, boot = boot),
  noa = med_iptw(filter(xs, p01xrkoa == 0), Y = hspss, A = v00edcv, M = smoker, C = age + agesq + p02sex + p02race, boot = boot),
  all = med_iptw(xs, Y = hspss, A = v00edcv, M = smoker, C = age + agesq + p02sex + p02race, boot = boot)
)

res_drk <- list(
  roa = med_iptw(filter(xs, p01xrkoa != 0), Y = hspss, A = v00edcv, M = drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  soa = med_iptw(filter(xs, p01sxkoa != 0), Y = hspss, A = v00edcv, M = drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  noa = med_iptw(filter(xs, p01xrkoa == 0), Y = hspss, A = v00edcv, M = drnkamt, C = age + agesq + p02sex + p02race, boot = boot),
  all = med_iptw(xs, Y = hspss, A = v00edcv, M = drnkamt, C = age + agesq + p02sex + p02race, boot = boot)
)