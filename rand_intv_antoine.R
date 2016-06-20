fmatch <- function(v, x, l = NULL, mat = FALSE){
# replaces syntax of type [cbind(1:nrow(dat),match(dat[[M]], levels(dat[[M]])))]
  if(!is.factor(x)) stop("Must be factor variable")
  if(is.null(l)) l <- length(x)
  
  if(length(levels(x)) == 2) v <- cbind(1-v, v)

  if(mat){return(cbind(1:l,match(x, levels(x))))
  }else return(v[cbind(1:l,match(x, levels(x)))])
}


rint_med.mkdata <- function(orig_dat, X, M, Y, C = "", L){
  ### replicate C with parentheses to allow arbitrary syntax in confounder specification (like equations or functions)
  if(C != "") C <- paste0("(",C,")")
  
  rint_med.mkform <- function(X, C, L, M){
    # Creates formulas to be passed to propensity score models
    a_form <- as.formula(paste0(X, "~", paste0(C, collapse = "+")))
    m_form <- as.formula(paste0(M, "~", paste(X, paste0(L, collapse = "+"), paste0(C, collapse = "+"), sep = "+")))
    l_form <- as.formula(paste0(L, "~", X, "+", paste0(C, collapse = "+")))
    return(list(a_form = a_form, m_form = m_form, l_form = l_form))
  }
  rint_med.mkmods <- function(dat, f){
    # Makes propensity score models from formulas in rint_med.mkform
    a_mod <- multinom(data = dat, f$a_form, trace = FALSE)
    m_mod <- multinom(data = dat, f$m_form, trace = FALSE)
    l_mod <- multinom(data = dat, f$l_form, trace = FALSE)
    return(list(a_mod = a_mod, m_mod = m_mod, l_mod = l_mod))
  }
  rint_med.mkden <- function(dat, mods, X, L){
    # Makes denominator 
    pmlac <- predict(object = mods$m_mod, type = "probs") %>% fmatch(dat[[M]])
    pac <- predict(object = mods$a_mod, type = "probs") %>% fmatch(dat[[X]])
    return(list(pmlac = pmlac, pac = pac, den = pac * pmlac))
  }
  rint_med.copy_data <- function(dat, X){
    ## One copy for each value in support of X
    new <- do.call(rbind, lapply(levels(dat[[X]]), function (i) mutate(dat, astar = i)))
    new$astar <- factor(new$astar, labels = levels(dat[[X]]))
    return(new)
  }
  
  rint_med.mknum <- function(dat, M, L, X, mmod, lmod){
    ## get p( l | a*,c) for each value in support of L
    plasc <- predict(object = lmod, newdata = mutate_(dat, .dots = setNames(list(~astar), X)), type = "probs")
    if(is.null(dim(plasc))) plasc <- cbind(1-plasc, plasc)
    
    ## get p( m | l,a*,c) for each value in support of L
    pmlasc <- lapply(levels(dat[[L]]), function (i) 
                      predict(object = mmod, newdata = mutate_(dat, .dots = setNames(list(~astar), X)) %>%
                                                       mutate_(.dots = setNames(list(~i), L)),
                              type = "probs") %>% fmatch(dat[[M]])
                     ) %>% do.call(cbind, .)
    
    ## sum p( m | l,a*,c ) * p( l | a*,c ) across levels of L
    num <- apply((plasc * pmlasc), 1, sum) 
    
    
    #Same pmlasc as AK
    pmlasc_ak <- predict(object = mmod, 
                       newdata = mutate_(dat, .dots = setNames(list(~astar), X)),
                       type = "probs") %>% fmatch(dat[[M]])
    num_ak <- apply((plasc * pmlasc_ak), 1, sum)
    return(list(plasc = plasc, pmlasc = pmlasc, num = num, num_ak = num_ak))
  }
  
  rint_med.mkweight <- function(dat){
    w <- dat$num/dat$den
    return(w)
  }
  #To get same weights as AK
  rint_med.mkweight_ak <- function(dat){
    w <- dat$num_ak/dat$den
    return(w)
  }
  
  ## create models
  flist <- rint_med.mkform(X, C, L, M)
  mlist <- rint_med.mkmods(orig_dat, flist)
  
  ## calculate denominator
  den_items <- rint_med.mkden(orig_dat, mlist, X, L)
  orig_dat$ipw_conf <- 1/den_items$pac
  orig_dat$den <- den_items$den
  
  ## duplicate dataset
  an_dat <- rint_med.copy_data(orig_dat, X)
  
  ## calculate numerator
  num_items <- rint_med.mknum(an_dat, M, L, X, mmod = mlist$m_mod, lmod = mlist$l_mod)
  an_dat$num <- num_items$num

  ## calculate weight (regular)
  an_dat$w <- rint_med.mkweight(an_dat)
  
  
  
  # get same weights as AK
  an_dat$num_ak <- num_items$num_ak
  an_dat$w_ak <- rint_med.mkweight_ak(an_dat)
    
  return(list(
    models = mlist,
    an_dat = an_dat,
    pmlac = den_items$pmlac,
    pac = den_items$pac,
    plasc = num_items$plasc,
    pmlasc = num_items$pmlasc,
    den = an_dat$den,
    num = an_dat$num,
    w = an_dat$w,
    
    #same weights as AK
    w_ak = an_dat$w_ak,
    num_ak = num_items$num_ak
  ))
}

rint_med.decompose <- function(dat, Y, X, astar = "astar"){
  ref <- levels(dat[[astar]])[1]
  m_te <- lm(data = dat, as.formula(paste0(Y, "~", X)), weights = ipw_conf)
  m_ter <- lm(data = dat[dat[[astar]] == dat[[X]],], as.formula(paste0(Y, "~", astar)), weights = w)
  m_nder <- lm(data = dat[dat[[astar]] == ref,], as.formula(paste0(Y, "~", X)), weights = w)
  m_nier_prollywrong <- lm(data = dat[dat[[X]] != ref,], as.formula(paste0(Y, "~", astar)), weights = w)
  nier <- lapply(levels(dat[[X]])[-1], function(i) lm(data = dat[dat[[X]] == i,], 
                                                     as.formula(paste0(Y, "~", astar)), weights = w)$coefficients[-1]
                ) %>% unlist
  
  return(list(
    models = list(m_te = m_te, m_ter = m_ter, m_nder = m_nder,
                  m_nier_prollywrong = m_nier_prollywrong),
    nder = m_nder$coefficients[-1],
    nier = nier,
    nier_prollywrong = m_nier_prollywrong$coefficients[-1],
    ter = m_ter$coefficients[-1],
    te = m_te$coefficients[-1]
  ))
}

rint_med2 <- function(dat, Y, X, astar = "astar"){
  ref <- levels(dat[[astar]])[1]
  
  return(unlist(nie))
}

  