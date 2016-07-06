source("rand_intv_alex.R")
source("rand_intv_antoine.R")

dgm <- function(samp_size){
  # simulate according to dag/sem given by
  # c <- {eps}
  # a <- {c, eps}
  # l <- {c, a, eps}
  # m <- {c, a, l, eps}
  # y <- {c, a, l, m, eps}
  c <- rbinom(samp_size, 1, 0.5)
  a <- rbinom(samp_size, 1, expit(-0.5 + c))
  l <- rbinom(samp_size, 1, expit(-0.5 + -0.5 + c + a))
  m <- rbinom(samp_size, 1, expit(-0.5 + -0.5 + -0.5 + c + a + l))
  y <- rnorm(samp_size, 0, 1)-0.5 + -0.5 + -0.5 + -2*0.5 + c + a + l + 2*m
  for(i in c("c", "a", "l", "m") ) assign(i, factor(get(i)))
  return(data.frame(id=1:samp_size, c=c, a=a, m=m, l=l, y=y))
}
orig_dat <- dgm(100000)

###### ALEX ################################################
an_dat <- mk.an.dat(orig_dat)

# randomized interventional analogue of natural direct/indirect effects
# NDE-R
m_nde <- glm(y ~ a_obs, weights = ipw, data=an_dat[an_dat$a==0,])

# NIE-R
m_nie <- glm(y ~ a, weights = ipw, data=an_dat[an_dat$a_obs==1,])

# randomized interventional analogue of the total effect 
m_ter = glm(y ~ a_obs, weights = ipw, data = an_dat[an_dat$a==an_dat$a_obs,])

# total effect (not the same)
m_te <- glm(y ~ a_obs, weights = ipw_conf, data = an_dat[an_dat$a==1,])
m_nde$coefficients[2]
m_te$coefficients[2] # TE
m_ter$coefficients[2] # TER
m_nde$coefficients[2] + m_nie$coefficients[2] #NDE-R + NIE-R



###### ANTOINE  ##############################################
orig_dat2 <- orig_dat

for(i in names(orig_dat2)) if(length(unique(orig_dat2[[i]])) < 3) orig_dat2[[i]] <- factor(orig_dat2[[i]])
rint_items <- rint_med.mkdata(orig_dat2, Y = "y", X = "a", C = "c", M = "m", L = "l")

#results following AB implementation
results <- rint_med.decompose(rint_items$an_dat, Y = "y", X = "a")

#results following AK implementation
results_ak <- rint_med.decompose(rint_items$an_dat %>% mutate(w = w_ak), Y = "y", X = "a")


##############################################################


rint_res.an <- rint_med.boot(filter(xs, p01sxkoa != 0), A = X, M = M, Y = Y, C = C, L = L,  boot = 100, alex = FALSE)

ssizes <- c(100, 1000, 10000, 100000, 300000, 500000, 1000000)

reslist <- lapply(ssizes, function(i){
    tmpd <- rint_med.mkdata(dgm(i), Y = "y", X = "a", C = "c", M = "m", L = "l")
    res <- rint_med.decompose(tmpd$an_dat, Y = "y", X = "a")
    c(res$ter, res$te, res$nier, res$nder)
  }) %>% do.call(rbind,. ) %>% unname

resdat <- reslist %>% cbind(ssizes) %>% data.frame()
colnames(resdat) <- c("ter", "te", "nier", "nder", "N")
resdat$sumr <- resdat$nder + resdat$nier




