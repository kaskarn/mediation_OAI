source("rand_intv_alex.R")
source("rand_intv_antoine.R")

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