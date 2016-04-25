#--------------------#
# Internal functions #
#--------------------#

# Internal functions that should not be exported

# Index

mleWrapper = function( dat, mle_fn, st_fn,
                       grad_fn = NULL, method = 'Nelder-Mead',
                       priors = NULL, SE = T, emStop = 20,
                       ... ) {
  # Purpose:
  # Forthcoming
  # Arguments:
  # Forthcoming
  # Returns:
  # Forthcoming

  # Generate starting values that don't produce invalid
  # log-likelihoods
  chk = -Inf
  em = 0
  while( chk == -Inf ) {
    startVal = st_fn()
    chk = try( mle_fn( startVal, dat, priors ), silent = T )
    if ( is.character( chk ) ) chk = -Inf
    em = em + 1
    if (em == emStop) stop("No feasible starting values were found",
                           call. = FALSE)
  }

  # Carry out optimization using function from 'stat' package
  results = try( optim( startVal, mle_fn,
                        dat = dat, priors = priors,
                        method = method, hessian = SE,
                        gr = grad_fn,
                        control = list( fnscale = -1, maxit = 5000, ... )
  ), silent = T )

  # Initialize output
  out = list(
    param = NA,
    logLik = NA,
    startVal = NA,
    convergenceCheck = NA,
    hessianMatrix = NA,
    optimOutput = NA
  )

  # If no errors, save output
  if ( is.list( results ) ) {
    out$param = results$par
    out$logLik = results$value
    out$startVal = startVal
    out$convergenceCheck = results$convergence
    if (SE) out$hessianMatrix = results$hessian
  }

  out$optimOutput = results

  return( out )
}


AICc_val = function( logLik, k, n ) {

  # Calculate information criterion
  aic = 2*k - 2*logLik
  # Apply correction for small samples
  # Note that as sample size grows, correction goes to 0
  aicc = aic + ( ( 2*k*( k+1 ) )/( n-k-1 ) )

  return( aicc )
}

BIC_val = function( logLik, k, n ) {

  # Calculate information criterion
  bic = -2*logLik + k*log(n)

  return( bic )
}

# Lookup - 06
covCreate = function(mat) {
  # Purpose:
  # Creates a single variable with a set of unique levels based on
  # a set of covariates, the number of all possible combinations for each
  # of the covariates
  # Arguments:
  # mat - A matrix of the covariates of interest
  # Returns:
  # A vector specifying the corresponding combination level for each
  # observation

  # Determine the different levels for each covariate
  lstLevel = lapply(as.data.frame(mat),unique)

  # Determine the possible combinations of the different covariates
  unq = expand.grid(lstLevel)
  Levels = 1:nrow(unq)

  # Output variable
  out = numeric( nrow(mat) )

  for ( k in Levels ) {

    for ( n in 1:nrow(mat) ) {
      if ( sum(mat[n,] == unq[k,])==ncol(mat) ) out[n] = k
    }
  }

  out
}

epanechnikovKernel = function(x) {

  k = 3/(4*sqrt(5))
  p1 = 1 - (x^2)/5
  p1[ abs(x) > sqrt(5) ] = 0

  k*p1
}

smoothing = function(x,num = .9) {

  n = length(x)
  s = sd(x)
  IQR = diff(quantile(x,prob=c(.25,.75)))

  out = ( num/(n^(1/5)) )*min(s,IQR/1.349)
  out
}

createIncrement = function(x) {
  # Purpose:
  # Forthcoming
  # Arguments:
  # x
  # Output:
  # Forthcoming

  # Determine the total number of original values
  curVal = sort( unique( x ) )
  newVal = 1:length(curVal) # Create regular increments
  new.x = numeric( length( x ) )
  for (i in 1:length(curVal)) {
    sel = x == curVal[i]
    new.x[sel] = newVal[i]
  }

  new.x
}
