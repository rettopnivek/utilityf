#--------------------#
# Internal functions #
#--------------------#

# Internal functions that should not be exported

# Index
# Lookup - 01:  AICc_val
# Lookup - 02:  BIC_val

# Lookup - 01
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

# Lookup - 02
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
