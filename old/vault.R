

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
