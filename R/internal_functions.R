#--------------------#
# Internal functions #
#--------------------#

# Internal functions that should not be exported

# Index
# Lookup - 01:  AIC_val
# Lookup - 02:  AICc_val
# Lookup - 03:  BIC_val
# Lookup - 04:  .ls.objects

# Lookup - 01
AIC_val = function( logLik, k, n ) {
  # Purpose:
  # Calculates Akaike's Information Criterion.
  # Arguments:
  # logLik - A summed log-likelihood value.
  # k      - The number of free model parameters.
  # n      - The number of observations.
  # Returns;
  # The AIC value.

  # Calculate information criterion
  aic = 2*k - 2*logLik

  return( aicc )
}

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

# Lookup - 04
# Include function from Dirk Eddelbuettel
# See http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
.ls.objects = function ( pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n = 5 ) {

  napply = function( names, fn)
    sapply( names, function(x) fn( get( x, pos = pos) ) )

  names = ls( pos = pos, pattern = pattern )
  obj.class = napply( names, function(x) as.character(class(x))[1] )
  obj.mode = napply( names, mode )
  obj.type = ifelse( is.na(obj.class), obj.mode, obj.class )
  obj.prettysize = napply( names, function(x) {
    format(utils::object.size(x), units = "auto")
  }
  )
  obj.size = napply( names, object.size )
  obj.dim = t( napply( names, function(x) as.numeric(dim(x))[1:2]) )
  vec = is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] = napply(names, length)[vec]
  out = data.frame( obj.type, obj.size, obj.prettysize, obj.dim )
  names(out) = c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out = out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out = head(out, n)
  out
}
