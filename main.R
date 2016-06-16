library(magrittr); library(mice)
library(foreign); library(dplyr)
library(reshape2); library(ipw)
library(foreach); library(nnet)

#load data and format
source("load_data.R")

#### utility functons
fmatch <- function(v, x, l = NULL, mat = FALSE){
  #replaces syntax of type [cbind(1:nrow(dat),match(dat[[M]], levels(dat[[M]])))]
    if(!is.factor(x)) stop("Must be factor variable")
    if(is.null(l)) l <- length(x)
    
    if(length(levels(x)) == 2) v <- cbind(1-v, v)
  
    if(mat){return(cbind(1:l,match(x, levels(x))))
    }else return(v[cbind(1:l,match(x, levels(x)))])
}
whichtimes <- function (stub, data = df, plen = 3){
  ## Summarizes data availability (for exploration)
  grep(paste0(stub, "$"), colnames(data), value = TRUE) %>%
    substring(., 2, 3) %>% as.numeric(.)
}

#### Functions based on vdw book
med.iptw <- function(dat, X, M, Y, C = "", L, xlen = 0, mlen = 0, mlvl = NULL, stable = FALSE){
  ### WORK IN PROGRESS
  if(xlen == 0) if(is.factor(dat[[X]])) xlen <- levels(dat[[X]]) %>% length - 1 else stop("X must be factor")
  if(mlen == 0) if(is.factor(dat[[M]])) mlen <- levels(dat[[M]]) %>% length - 1 else stop("M must be factor")
  
  #replicate input with parentheses
  for(i in c("X","M","Y","C", "L")) assign(paste0(i, "2"), paste0("(",get(i),")"))
  
  N <- nrow(dat)
  #p(a)
  pa <- (table(dat[[X]])/N)[match(dat[[X]], levels(dat[[X]]))]
  pa_c <- multinom(data = dat, as.formula(paste0(X, "~", paste0(C2, collapse = "+"))), trace = FALSE) %>% 
          predict(., newdata = dat, "probs") %>% fmatch(dat[[X]])
  pm_a <- multinom(data= dat, as.formula(paste0(M, "~", X)), trace = FALSE) %>% 
          predict(., newdata = dat, "probs") %>% fmatch(dat[[M]])
  pm_acl <- multinom(data= dat, as.formula(paste0(M, "~", paste(X, paste0(C2, collapse = "+"), paste0(L2, collapse = "+"), sep = "+"))), trace = FALSE) %>% 
            predict(., newdata = dat, "probs") %>% fmatch(dat[[M]])
  w <- pa*pm_a/pa_c/pm_acl
  
  #run MSM
  msm <- lm(data = dat, get(Y) ~ get(X) + get(M) + get(X)*get(M), weights = w)
  
  #separate output
  num <- cumsum(c(xlen, mlen, xlen*mlen))
  
  #treat mlvl
  if(is.null(mlvl)) mlvl <- table(dat[[M]]) / N
  cde <- lapply(mlvl)
}
twoway.cont.cont <- function(df, X, M, Y, C = "", xlen = 0, mlen = 0, noint = FALSE, mlvl = NULL, delta = FALSE){
  ## Twoway decomposition indirect and direct effects  (VDW)
  #for continuous mediator and outcome

  
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
rint_med <- function(dat, X, M, Y, C = "", L, xlen = 0, mlen = 0, stable = FALSE, debug = FALSE){
  # Random interventional analogue to natural direct and indirect effects
  
  ###Create shorthand for base level of X
  ref <- levels(dat[[X]])[1]

  ###replicate input with parentheses. 
  #This allow to use arbitrary syntax in confounder specification (like equations or functions)
  #functionally, X, M, Y, C and L are the same as X2, Y2, M2, C2, L2 in most cases
  for(i in c("X","M","Y","C", "L")) assign(paste0(i, "2"), paste0("(",get(i),")"))
  
  ### Compute weights
  #[ ]_form: formula fed to mode
  #[ ]_mod: model of a/m/l, 
  #[ ]_pro: propensity score (conditional probability)
  #note 1: dplyr functions (mutate_) probably have faster datatable analogues, but I prefer dplyr grammar
  #note 2: requires fmatch() function defined above, which returns conditional probability of observed from matrix of 
  # predicted probabilities for all possible levels
  
  ## denominator
  # p(a | c)
  a_form <- as.formula(paste0(X2, "~", paste0(C2, collapse = "+")))
  if(debug) print(a_form)
  a_mod <- multinom(data = dat, a_form, trace = FALSE)
  a_pro <- predict(object = a_mod, newdata = dat, "probs") %>% fmatch(dat[[X]])
  
  # p(m |l, a, c)
  m_form <- as.formula(paste0(M2, "~", paste(X, paste0(L2, collapse = "+"), paste0(C2, collapse = "+"), sep = "+")))
  if(debug) print(m_form)
  m_mod <- multinom(data = dat, m_form1, trace = FALSE)
  m_pro <- predict(object = m_mod1, newdata = dat, "probs") %>% fmatch(dat[[M]])

  den <- a_pro * m_pro
  
  ## numerator
  # p(l | a*, c)
  l_form <- as.formula(paste0(L2, "~", X, "+", paste0(C2, collapse = "+")))
  if(debug) print(l_form)
  l_mod <- multinom(data = dat, l_form, trace = FALSE)
  l_pro <- predict(object = l_mod, newdata = mutate_(dat, .dots = setNames(list(~ref), X)), "probs")
  
  # sum over support of L: P(m | a*, c, l)P(l | a*, c)
  num <- lapply(levels(dat[[L]]), 
                  function (i) {
                      predict(
                        object = m_mod, #model of m from before
                        newdata = mutate_(dat, .dots = setNames(list(~i), L)) %>% #nonstandard evaluation to set L to current level in sum
                                  mutate_(., .dots = setNames(list(~ref), X)), #set A to unexposed (can be taken out of loop)
                        "probs"
                      ) %>% fmatch(dat[[M]]) * l_pro[,match(i, levels(dat[[L]]))] #P(m | a*, c, l)P(l | a*, c)
                  }
                ) %>% 
          do.call(cbind, .) %>% #forms a matrix of P(m | a*, c, l)P(l | a*, c) for each l in support of L 
          apply(., 1, sum) #sums across rows for marginal p
  
  #get weight, stabilize if needed (shouldn't not required)
  if(stable){ 
    dat$w <- num/den * ((table(dat[[X]])/nrow(dat))[match(dat[[X]], levels(dat[[X]]))])
  }else dat$w <- num/den  
  if(debug) mean(w, na.rm = TRUE) %>% print
  
  ### Direct effect
  nde <- lm(data = dat, as.formula(paste0(Y2, "~", X2)), weights = w)$coefficients[-1] #we don't save the intercept coeff
  
  ### Indirect effect
  #duplicate dataset with 0/1 indicator
  dat <- rbind(mutate(dat, astar = 0), mutate(dat, astar = 1)) %>% mutate(astar = factor(astar))
  #get NIEs in models for each nonrefent level of X
  nie <- lapply(levels(dat[[X]])[-1], function (i) lm(data = dat[dat[[X]] == i,], 
                                    as.formula(paste0(Y2, "~ astar")), 
                                    weights = w)$coefficients[2]) #save astar coeff
  
  ### output results
  nie <- unlist(nie)
  names(nie) <- names(nde) <- levels(dat[[X]])[-1]
  return(list(mod_exp = a_mod,
              mod_med = m_mod,
              nde = nde, 
              nie = nie))
}

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
test <- rint_w(xs, X = "v00edcv", M = "obcat", Y = "hspss", C = confounders, L = "smoker")

#######################################

# 
# ### tests of regresson methdos
# 
# #test simple regression-based methods
# X <- "bmi"
# M <- "hsmss"
# C <- confounders
# 
# test <- twoway.cont.cont(xs, X = "exposure", M = "bmi", Y = "hspss", C = confounders, mlvl = 25:30, noint = FALSE)
# test <- twoway.cont.cont(xs, X = "exposure", M = "bmi", Y = "hspss", C = confounders, mlvl = 25:30, noint = FALSE)
# ########################################   