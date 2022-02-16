
PReMo_model <- function(parfree, parset, schedule) {
  
  # we want all the parameters in 1 vector:
  p <- c(parfree, parset)
  
  # so that we can read them out as variables here:
  var_u     <- p['var_u']
  var_p     <- p['var_p']
  var_v     <- p['var_v']
  eta_p     <- p['eta_p']
  eta_v     <- p['eta_v']
  beta_psat <- p['beta_psat']
  beta_vsat <- p['beta_vsat']
  K         <- p['K']
  
  # initialize the state:
  x_pt    <- 0
  x_vt    <- 0
  beta_pt <- 0
  beta_vt <- 0
  Gt      <- 0 # this will be kept 0, could increase for explicit re-aiming
  
  # empty vectors:
  reachdev  <- c()
  proprecal <- c()

  for (trial in c(1:dim(schedule)[1])) {
    
    # Marius' added line:
    x_vt <- switch(schedule$trialtype[trial],
                   x_pt + schedule$rotation[trial],
                   schedule$rotation[trial])
    
    # Eq. 1
    # integrated proprioceptive position:
    xI_pt <- ( (var_u/(var_u+var_p)) * x_pt ) + ( (var_p/(var_u+var_p)) * Gt )
    
    # Eq. 2
    # integrated visual position:
    xI_vt <-  ( (var_u/(var_u+var_v)) * x_vt ) + ( (var_p/(var_u+var_v)) * Gt )
    
    # Eq. 3
    # perceived proprioceptive position:
    xper_pt <- xI_pt + beta_pt
    
    # Eq. 4
    # perceived visual position:
    xper_vt <- xI_vt + beta_vt
    
    # Eq. 5
    # reported proprioceptive shift:
    beta_pt <- min(eta_p*(xI_vt - xI_pt), beta_psat)
    
    # Eq. 6
    # reported visual shift:
    beta_vt <- min(eta_p*(xI_vt - xI_pt), beta_vsat)
    
    
    # store stuff in vectors:
    reachdev  <- c(reachdev,  x_pt)
    proprecal <- c(proprecal, beta_pt)
    
    # Eq. 8 (and 7)
    # reach direction on next trial:
    x_pt <- x_pt + (K*(Gt - xper_pt))
    
  }
  
  return(data.frame(reachdev, proprecal))
  
}



PReMo_errors <- function(parfree, parset, schedule, dataset) {
  
  # this will compare model output to reaches, or proprioception, or both
  # but not neither
  
  if ( !any(c('reachdev', 'proprecal') %in% names(dataset)) ) {
    cat('can not evaluate model errors if there is no data\n')
  }
  # can do more checks, such as that if these variables are not null they are: 
  # - some kind of numeric vector
  # - with length equal to schedule
  
  # run model with current parameters:
  prediction <- PReMo_model(parfree, parset, schedule)
  
  # collect errors here:
  errors <- c()
  
  for (varname in c('reachdev', 'proprecal')) {
    if (varname %in% names(dataset)) {
      errors <- c(errors, prediction[,varname] - dataset[,varname])
    }
  }
  
  # remove NAs
  errors <- errors[!is.na(errors)]
  
  return(sum(errors^2))
  
}


library(optimx)

PReMo_fit <- function(dataschedule, 
                      parset=c(),
                      npoints=6,
                      maxgridsize=250,
                      gridtop=5,
                      parlimits = NULL) {
  
  # split dataschedule:
  dataset  <- dataschedule[,which(names(dataschedule) %in% c('reachdev','proprecal'))]
  schedule <- dataschedule[,which(names(dataschedule) %in% c('rotation','trialtype'))]
  
  
  # start search grid for free parameters:
  searchseqs <- list()
  lower <- c()
  upper <- c()
  
  # if not set, these are the default limits:
  if (is.null(parlimits)) {
    parlimits <- data.frame('parameter'=c('K',  'eta_p', 'eta_v', 'var_u', 'var_p', 'var_v', 'beta_vsat', 'beta_psat'),
                            'lo'       =c(0.001, 0.001,  0.001,    0.001,    0.001,   0.001,   0.001,       0.001),
                            'hi'       =c(0.999, 0.999,  0.999,   20,       20,       5,      10,          15))
  }
  
  for (rown in c(1:dim(parlimits)[1])) {
    if (parlimits$parameter[rown] %in% names(parset)) { next() }
    lo <- parlimits[rown,'lo']
    hi <- parlimits[rown,'hi']
    if (lo == hi) {
      thisparval <- c(lo)
      names(thisparval) <- c(parlimits$parameter[rown])
      parset <- c(parset, thisparval)
    } else {
      step <- diff(c(lo,hi))/npoints
      searchseqs[[parlimits$parameter[rown]]] <- seq(lo+(step/2),hi-(step/2),step)
      lower <- c(lower, lo)
      upper <- c(upper, hi)
    }
  }
  
  print(parset)
  
  searchgrid <- expand.grid(searchseqs)
  
  if (dim(searchgrid)[1] > maxgridsize) {
    searchgrid <- searchgrid[sample(c(1:dim(searchgrid)[1]),size=maxgridsize),]
  }
  
  MSEs <- apply(searchgrid,FUN=PReMo_errors,MARGIN=c(1),parset=parset,schedule=schedule,dataset=dataset)
  
  topgrid <- searchgrid[order(MSEs)[c(1:gridtop)],]
  
  # do fitting on the best results from the search grid:
  topFits <- do.call("rbind",
                     apply( topgrid,
                            MARGIN=c(1),
                            FUN=optimx::optimx,
                            fn=PReMo_errors,
                            method='L-BFGS-B',
                            lower=lower,
                            upper=upper,
                            schedule=schedule,
                            parset=parset,
                            dataset=dataset) )
  
  winFit <- topFits[order(topFits$value)[1],]
  
  idx <- which(names(winFit) == 'value') - 1
  
  winpar <- as.numeric(winFit[c(1:(idx))])
  names(winpar) <- names(winFit[c(1:(idx))])
  
  return(c(winpar,parset))
  
}

PReMo_fastfit <- function(dataschedule, 
                          settings) {
  
  # split dataschedule:
  dataset  <- dataschedule[,which(names(dataschedule) %in% c('reachdev','proprecal'))]
  schedule <- dataschedule[,which(names(dataschedule) %in% c('rotation','trialtype'))]
  
  # if not set, these are the default limits:
  parlimits <- data.frame('parameter'=c('K',  'eta_p', 'eta_v', 'var_u', 'var_p', 'var_v', 'beta_vsat', 'beta_psat'),
                          'lo'       =c(0.001, 0.001,  0.001,    0.001,    0.001,   0.001,   0.001,       0.001    ),
                          'hi'       =c(0.999, 0.999,  0.999,   20,       20,       5,      10,          15        ))
  
  lower   <- c()
  upper   <- c()
  parset  <- c()
  parfree <- c()
  for (parname in parlimits$parameter) {
    pl_idx <- which(parlimits$parameter == parname)
    st_idx <- which(settings$parameter == parname)
    
    if (settings$fixed[st_idx] == TRUE) {
      parset[parname] <- settings$value[st_idx]
    } else {
      thispar <- c()
      thispar[parname] <- settings$value[st_idx]
      parfree <- c(parfree, thispar)
      lower <- c(lower, parlimits[pl_idx,'lo'])
      upper <- c(upper, parlimits[pl_idx,'hi'])
    }
    
  }
  
  fit <- optimx::optimx(par=parfree,
                        fn=PReMo_errors,
                        method='L-BFGS-B',
                        lower=lower,
                        upper=upper,
                        schedule=schedule,
                        parset=parset,
                        dataset=dataset)
  
  return(fit)
  
}
