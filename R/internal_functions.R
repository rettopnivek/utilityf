#--------------------#
# Internal functions #
#--------------------#

# Internal functions that should not be exported

# Index
# Lookup - 01:  mleWrapper
# Lookup - 02:  AICc_val
# Lookup - 03:  BIC_val

# Lookup - 01
mleWrapper = function( dat, mle_fn, st_fn,
                       grad_fn = NULL, method = 'Nelder-Mead',
                       priors = NULL, SE = T, emStop = 20,
                       ... ) {
  # Purpose:
  # An internal function that uses R's base 'optim' function
  # to carry out maximum likelihood optimization.
  # Arguments:
  # dat     - an R object with the data to be fitted. Passed on to
  #           the 'mle_fn' function.
  # mle_fn  - a R function that must take three variables as inputs:
  #           1) a vector of parameters;
  #           2) the data (i.e. 'dat');
  #           3) a set of prior values (can be NULL).
  # st_fn   - a R function that generates a set of starting values for
  #           the parameter vector (must take no inputs).
  # grad_fn - an optional R function that calculates a vector of the
  #           first derivatives for the parameters.
  # method  - the method the 'optim' function should use
  # priors  - the prior values to pass on to 'mle_fn'.
  # SE      - a logical value indicating if standard errors should be
  #           calculated.
  # emStop  - the maximum number of times generation of starting values
  #           is attempted.
  # ...     - Addtional variables to pass on to the control list for
  #           'optim'
  # Returns:
  # A list consisting of...
  # 1) param            - the vector of best-fitting parameters
  # 2) logLik           - the maximum log-likelihood value
  # 3) startVal         - the set of starting values
  # 4) convergenceCheck - the optim code for convergence
  # 5) hessianMatrix    - if SE == T, the hessian matrix that was calculated
  # 6) optimOutput      - the list that 'optim' outputs

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
                        control = list( fnscale = -1, ... )
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
# Example of function:
# mle_fn = function( par, dat, priors ) {
#   logLik = sum( dnorm( dat, par[1], exp( par[2] ), log = T ) );
#   prs = dnorm( par[1], priors[1,1], priors[1,2], log = T );
#   prs = prs + dgamma( exp( par[2] ), priors[2,1], priors[2,2], log = T );
#   return( logLik + prs ) }
# st_fn = function() { out = runif( 2, c(-3,.1), c(3,3) ); out[2] = log( out[2] ); out }
# dat = rnorm( 100, 10, 3 ); priors = rbind( c(0,4), c(2,2) );
# tst = mleWrapper( dat, mle_fn, st_fn, priors = priors )

# Lookup - 02
AICc_val = function( logLik, k, n ) {
  # Purpose:
  # Calculates Akaike's Information Criterion with a correction
  # for small sample sizes.
  # Arguments:
  # logLik - A summed log-likelihood value.
  # k      - The number of free model parameters.
  # n      - The number of observations.
  # Returns;
  # The AICc value.

  # Calculate information criterion
  aic = 2*k - 2*logLik
  # Apply correction for small samples
  # Note that as sample size grows, correction goes to 0
  aicc = aic + ( ( 2*k*( k+1 ) )/( n-k-1 ) )

  return( aicc )
}

# Lookup - 03
BIC_val = function( logLik, k, n ) {
  # Purpose:
  # Calculates the Bayesian Information Criterion.
  # Arguments:
  # logLik - A summed log-likelihood value.
  # k      - The number of free model parameters.
  # n      - The number of observations.
  # Returns;
  # The BIC value.

  # Calculate information criterion
  bic = -2*logLik + k*log(n)

  return( bic )
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
