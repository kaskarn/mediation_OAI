# Alex Keil
# rand_intv_effectdecomp.R
# 6/16/2016
# quick simulation example of 
# randomized interventional analogue of natural direct/indirect effects
# see also: 
# [1] T. J. Vanderweele, S. Vansteelandt, and J. M. Robins. 
#  Effect decomposition in the presence of an exposure-induced mediator-outcome confounder. 
#  Epidemiology, 25(2):300–6, Mar 2014.
# [2] T. VanderWeele. Explanation in causal inference: methods for mediation and interaction. 
#  Oxford University Press, 2015.


#helpers
expit <- function(mu) 1/(1+exp(-mu))

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
  return(data.frame(id=1:samp_size, c=c, a=a, m=m, l=l, y=y))
}

mk.ipw.mods <- function(dat){
  # make IPW based on approach 3 of Van. & Van 2014
  amod <- glm(a ~ c, data=dat, family=binomial)
  mmod <- glm(m ~ l + a + c, data=dat, family=binomial)
  lmod <- glm(l ~ a + c, data=dat, family=binomial)
  return(list(amod=amod, mmod=mmod, lmod=lmod))
}

get.denoms <- function(modlist, dat){
  # get denominator components of IPW
  # approach 3 of Van. & Van 2014
  # p(a|c)
  ps <- predict(modlist$amod, type='response')
  pac <- ifelse(dat$a, ps, 1-ps) #### 
  # p(m|lac)
  pm <- predict(modlist$mmod, type='response')
  pmlac <- ifelse(dat$m, pm, 1-pm)
  return(list(pac=pac, pmlac=pmlac))
}


copy.data <- function(dat, a_levels=2){
  # make copies of dataset with pseudo exposure
  dlist <- c()
    dat$a_obs <- dat$a
  for(setA in 0:1){
    dat$a <- setA
    dlist <- c(dlist, paste0('dat', setA))
    assign(paste0('dat', setA), dat)
  }
  #########################
  # todo generalize this to more exposure levels
  rbind(dat0, dat1)
}


mk.an.dat <- function(dat){
  #make analytic dataset
  #create ipws, output dataset for analysis
  #p(m|l,a*,c)
  i_mods <- mk.ipw.mods(dat)
  d_list <- get.denoms(i_mods, dat)
  dat$pac <- d_list$pac
  #########################
  # todo generalize this to more exposure levels
  dat$ipw_conf <- 1/dat$pac                         ## ////!\\\\ changed to inverse of pac, =/= ifelse(a, pac, 1-pac)
  #########################
  dat$pmlac <- d_list$pmlac
  # make copies
  an_dat = copy.data(dat)
  #p(m|l,a*,c)
  pmnum = predict(i_mods$mmod, newdata=an_dat, type='response')
  pmlasc = ifelse(an_dat$m, pmnum, 1-pmnum)
  #p(l|a*,c)
  plasc = predict(i_mods$lmod, newdata=an_dat, type='response')
  # sum the terms in the numerator
  num = pmlasc*plasc + pmlasc*(1-plasc)             ### ////!\\\\ factoring p(m | la*c) out of the sum over support of L, but conditional on L. Check?
  an_dat$ipw = num/(an_dat$pac*an_dat$pmlac)
  an_dat
}