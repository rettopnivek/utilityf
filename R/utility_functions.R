#----------------------------#
# Assorted utility functions #
#----------------------------#

# Package development
# library(devtools)
# library(roxygen2)

# Index
# Lookup - 01:  Logit Function
# Lookup - 02:  Logistic Function
# Lookup - 03:  Plotting Range Function
# Lookup - 04:  Wrapper for Maximum Likelihood Estimation
# Lookup - 05:  Function to Draw Ellipses
# Lookup - 06:  Error Function
# Lookup - 07:  Highest Density Interval Estimator
# Lookup - 08:  Standard Error of the Mean
# Lookup - 09:  Generate a Blank Plot
# Lookup - 10:  Calculate Standard Errors and Confidence Intervals
# Lookup - 11:  Conversion between Degrees and Radians
# Lookup - 12:  Conversion to Cartesian Coordinates
# Lookup - 13:  Create a Design Matrix
# Lookup - 14:  Softmax Function
# Lookup - 15:  Reverse of the Softmax Function

# Lookup - 01
#' Logit Function
#'
#' This function calculates the logit (log of the odds) of a set of probabilities.
#'
#' @param p a set of probabilities (0 <= p <= 1).
#' @return Returns a set of values that now lie between -Inf and Inf.
#' @export
#' @examples
#' logit( c(.5,.1,.9,0,1) )

logit = function(p) {
  p[p<0 | p>1] = NA
  log(p/(1-p))
}

# Lookup - 02
#' Logistic Function
#'
#' This function applies the logistic function to a set of values
#'
#' @param x a set of values ( -Inf <= x <= Inf).
#' @return Returns a set of probabilities.
#' @keywords logistic
#' @export
#' @examples
#' logistic( c(0,-2.197225,2.1972255,-Inf,Inf) )

logistic = function(x) {
  1/(1+exp(-x))
}

# Lookup - 03
#' Plotting Range Function
#'
#' This function determines the lower and upper boundaries for
#' a vector of values given a specified sub-interval for
#' plotting purposes.
#'
#' @param int the desired sub-interval over the range.
#' @param dat the vector of values over which to determine the boundaries.
#' @return Returns a lower and upper boundary, multiples of the
#' specified sub-interval.
#' @export
#' @examples
#' d = density( rnorm(100) )
#' x.limit = lowerUpper( .5, d$x )
#' y.limit = lowerUpper( .1, d$y )
#' plot( x.limit, y.limit, type='n', xlab='Values', ylab='Density' )
#' lines(d$x,d$y)

lowerUpper = function(int,dat) {

  ll = int*floor(min(dat)/int)
  ll = c(ll,int*ceiling(max(dat)/int))

  ll
}

# Lookup - 04
#' Wrapper for Maximum Likelihood Estimation
#'
#' This function allows multiple runs of the optim function
#' with random starting values to better control for local
#' maxima/minima.
#'
#' @param dat the data to be fitted (formats can vary).
#' @param st.fn a function to generate a set of starting values
#'   (must take no parameters and output must match the input length
#'   for the mle.fn function).
#' @param mle.fn a function to calculate the negative of the sum of the log-likelihood
#'   (must take two inputs: 1. a vector of parameters, 2. the data).
#'   Note that other estimation methods (e.g. least-squares) can be used
#'   instead if the function mle.fn is defined appropriately.
#' @param nRep the number of times to run the estimation routine.
#'   Running the routine multiple times with dispersed starting values
#'   increases the chances of finding the true maximum/minimum instead of
#'   local maxima/minima.
#' @param unDef the value returned by mle.fn indicating undefined values
#'   (default is Inf).
#' @param ... additional parameters for the optim function. For example
#'   additional options can be passed into the parameter 'control' as a list:
#'   control = list( maxit = 1000, fnscale = -1 ), which increases the
#'   number of iterations and use maximization instead of minimization.
#' @return Returns a list with 4 elements:
#'   MLE - a list of the successful outputs from the optim function;
#'   start.val - a matrix of the starting values for each successful
#'     run;
#'   logLikSum - a vector of the mle.fn function value for each
#'     successful run;
#'   run.time - The amount of time the algorithm took.
#' @export
#' @examples
#' mle.fn = function(par,dat) -sum(dnorm(dat,par[1],exp(par[2]),log=T)) # Likelihood function
#' st.fn = function() runif( 2, c(50,0),c(150,3) ) # Starting values function
#' dat = rnorm( 100, 100, 15 ) # Generate data
#' results = MLE( dat, st.fn, mle.fn )
#' sel = which( results$logLikSum==min(results$logLikSum)) # Select the minimum
#' results$MLE[[sel]] # The results
#' par = results$MLE[[sel]]$par; round( c(par[1],exp(par[2])) ) # The paramaters

MLE = function(dat,st.fn,mle.fn,nRep=4,unDef=Inf,...) {

  # Check that 'stats' package is installed
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("The 'stats' package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Track running time
  st.time = Sys.time()

  # Create set of variables to store results
  logLikSum = rep(NA,nRep)
  start.val = st.fn() # Determine initial length
  start.val = matrix(NA,nRep,length(start.val))
  results = c()
  for (i in 1:nRep) results=c(results,list(NULL))

  # Loop through repetitions
  inc = 1
  for (nr in 1:nRep) {

    # Generate random starting values that
    # don't produce invalid likelihoods
    chk = unDef
    while (chk==unDef) {
      st = st.fn()
      chk = mle.fn(st,dat)
    }

    # Set default value in case estimation fails
    res = NULL

    # Run estimation routine within a 'try' statement
    # so that function won't stop if there's an error
    res = try( optim(st,mle.fn,dat=dat,...) )
    # print(res)

    # Save results if no error occurred
    if (is.list(res)) {
      results[[inc]] = res
      start.val[inc,]=st
      logLikSum[inc]=res[['value']]
      inc = inc + 1
    }

  }

  # Remove any missing values
  chk = which(is.na(logLikSum))
  if (length(chk)>0) {
    results = results[-chk]
    start.val = start.val[-chk,]
    logLikSum = logLikSum[-chk]
  }

  run.time = Sys.time() - st.time

  # Output results as a list
  return(
    list(
      MLE = results,
      start.val = start.val,
      logLikSum = logLikSum,
      run.time = run.time
    )
  )
}

# Lookup - 05
#' Function to Draw Ellipses
#'
#' This function will draw an ellipse on a plot that
#' already exists.
#'
#' @param a the length of the x-axis vertice of the ellipse.
#' @param b the length of the y-axis vertice of the ellipse.
#' @param k the angle between the x-axis in the major vertice
#'   for the ellipse (default is 0).
#' @param Xc the center on the x-axis for the ellipse
#'   (default is 0).
#' @param Yc the center on the y-axis for the ellipse
#'   (default is 0).
#' @param deg if TRUE, assumes k is in degrees; otherwise,
#'   k is assumed to be in radians (default is TRUE).
#' @param lngth the number of points to use when plotting
#'   the ellipse (default is 100).
#' @param draw if TRUE, draws a polygon in the shape of the
#'   ellipse on an existing plot (default is TRUE).
#' @param ret if TRUE, returns a matrix with the x and y
#'   values for the ellipse (default is FALSE).
#' @param ... additional parameters for the polygon function.
#' @export
#' @examples
#' plot( c(-1,1), c(-1,1), type='n', xlab='X-axis', ylab='Y-axis' )
#' drawEllipse( .2, .2 )
#' ex = drawEllipse( .5, .2, Xc = -.5, Yc = .5, k = 20, ret = T, col = 'green')
#' drawEllipse( .05, .3, Xc = .5, Yc = -.5, k = 290, col = 'red')

drawEllipse = function(a, b, k = 0, Xc = 0, Yc = 0,
                       deg = T, lngth = 100, draw = T,
                       ret = F,...) {

  if (deg) k = k * (pi/180) # Convert to radians

  t = seq(0,2*pi,length=lngth)

  x = Xc + a * cos(t) * cos(k) - b * sin(t) * sin(k)
  y = Yc + b * sin(t) * cos(k) + a * cos(t) * sin(k)
  polygon(x,y,...)

  if (ret) return( cbind(x,y) )
}

# Lookup - 06
#' Error Function
#'
#' Calculates the error function.
#'
#' @param x a vector of values on the real number line.
#' @return a vector of transformed values based on the error function.
#' @export
#' @examples
#' x11(); plot( c(-1,1), c(-3,3), type='n', xlab='x', ylab='erf(x)' )
#' x = seq(-3,3,length=100)
#' lines(x,erf(x))

erf = function(x) {
  2*pnorm( x*sqrt(2), 0, 1 ) - 1
}

# Lookup - 07
#' Highest Density Interval Estimator
#'
#' Estimates the highest density interval for a specified coverage
#' for a sample drawn from a distribution of interest. Useful for
#' determining credible intervals given a MCMC sample.
#'
#' @param sampleVec a vector of representative values from a probability
#'   distribution.
#' @param credMass a scalar between 0 and 1, indicating the mass
#'   within the credible interval that is to be estimated.
#' @return A vector containing the limits of the HDI.
#' @section References:
#' Kruschke, J. (2010). Doing Bayesian data analysis: A tutorial
#'   introduction with R. Academic Press.
#' @export
#' @examples
#' round( qbeta( c(.025,.975), 10, 10 ), 3 ) # True values
#' round( hdi( rbeta(1000, 10, 10 ) ), 3 ) # Should be close

hdi = function( sampleVec, credMass=0.95 ) {
  sortedPts = sort( sampleVec )
  ciIdxInc = floor( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

# Lookup - 08
#' Standard Error of the Mean
#'
#' This function calculates the standard error of the mean.
#'
#' @param x a vector of values.
#' @return Returns the estimated standard error of the mean
#'   based on the values of x
#' @export
#' @examples
#' 1/sqrt(1000) # True value
#' sem( rnorm(1000) )

sem = function(x) {
  sd(x)/sqrt( length(x) )
}

# Lookup - 09
#' Generate a Blank Plot
#'
#' This function generates a completely blank plot.
#'
#' @param xDim the lower and upper boundaries for the
#'   x-axis (default values are [0,1] ).
#' @param yDim the lower and upper boundaries for the
#'   y-axis (default values are [0,1] ).
#'
#' @return A empty plot.
#' @export

blankPlot = function( xDim = c(0,1), yDim = c(0,1) ) {

  plot( xDim, yDim, type = 'n', ylab = ' ', xlab = ' ',
        xaxt = 'n', yaxt = 'n', bty = 'n' )

}

# Lookup - 10
#' Calculate Standard Errors and Confidence Intervals
#'
#' This function calculates the standard errors and confidence
#' intervals for a set of parameters extracted from an optim
#' object.
#'
#' @param fit an optim object in which the option Hessian == T.
#' @param alpha the range for the confidence interval (default
#'   is 0.95).
#'
#' @return A list consisting of two objects:
#'   SE - a vector with the standard errors for a set of parameters.
#'   CI - a matrix with two rows, the top row contains the lower
#'        intervals, the bottom row with the upper intervals.
#' @examples
#' mle.f = function( par, x ) { -sum( dnorm(x,par[1],par[2],log=T) ) }
#' x = rnorm(100,5,2)
#' res = optim( c(0,1), mle.f, x = x, hessian = T )
#' mle.SE( res )
#'
#' mle.f = function( par, x ) { sum( dnorm(x,par[1],par[2],log=T) ) }
#' x = rnorm(100,5,2)
#' res = optim( c(0,1), mle.f, x = x, hessian = T, control = list( fnscale = -1 ) )
#' mle.SE( res, alpha = .99, fnscale = -1 )
#' @export

mleSE = function( fit, alpha = .95, fnscale = 1 ) {

  fisher_info = solve(fnscale*fit$hessian)
  prop_sigma = sqrt(diag(fisher_info))
  prop_sigma = prop_sigma

  ub = qnorm( (1 - alpha)/2 )
  lb = qnorm( (1 - alpha)/2, lower.tail = F )

  upper<-fit$par+lb*prop_sigma
  lower<-fit$par+ub*prop_sigma

  return( list( SE = prop_sigma,
                CI = matrix(c(lower,upper),2,
                            length(prop_sigma),byrow=T) ) )
}


# Lookup - 11
#' Conversion between Degrees and Radians
#'
#' A function to convert degrees to radians and vice versa.
#'
#' @param degrees radians equal degrees(pi/180).
#' @param radians degrees equal 180*radians/pi.
#'
#' @return The converted values.
#' @examples
#' r = degreesRadians( degrees = 45 )
#' print( r )
#' d = degreesRadians( radians = .7853982 )
#' print( d )
#' @export

degreesRadians = function( degrees = NULL, radians = NULL ) {

  out = NULL

  if ( length( degrees )>0) {
    out = degrees*(pi/180)
  }

  if ( length( radians )>0) {
    out = radians/(pi/180)
  }

  out
}

# Lookup - 12
#' Conversion to Cartesian Coordinates
#'
#' A function to use the magnitude and angle of a 2D vector
#' and convert them to cartesian coordinates (assumes the
#' vector starts at [0,0] ).
#'
#' @param H a magnitude.
#' @param A an angle.
#' @param degrees if true, the angle is assumed to be in
#'   degrees.
#'
#' @return A set of cartesian coordinates.
#' @examples
#' cart = convertMagnitudeAngle( 2, 45 )
#' print( cart )
#' @export

convertMagnitudeAngle = function( H, A, degrees = T) {

  # Convert to radians
  if (degrees) A = degreesRadians( degrees = A )

  a = cos(A)*H; b = sin(A)*H

  c(a,b)
}

# Lookup - 13
#' Create Design Matrix
#'
#' A function that translates a vector of categories into different
#' types of design matrices.
#'
#' @param X a vector of categories
#' @param Levels the unique values of X. If left unspecified, the
#'   function attemps to extract, but the mappings will be based on order
#'   of appearance.
#' @param Mapping For 'Dummy' and "Effects' options, provides an
#'   additional way to indicates which unique values of X should be
#'   mapped to which columns of the design matrix.
#'   For the 'Coef' option, provides the corresponding weight values for
#'   the unique values of X.
#' @param type Current options include...
#'   \itemize{
#'     \item 'Dummy' -> Given K levels and N trials, returns a N x (K-1)
#'     design matrix using dummy coding (i.e. 0 and 1, with one variable).
#'     Useful for simple effects.
#'     \item 'Effects' -> Given K levels and N trials, returns a N x (K-1)
#'     design matrix using effects coding (i.e. -1 and 1, with one variable
#'     denoted solely by -1). Useful for comparisons against the grand mean.
#'     \item 'Intercept' -> Given K levels and N trials, returns a N x K
#'     design matrix in which each unique level has its own column (i.e. K
#'     unique intercepts).
#'     \item 'Coef' -> Given K levels and N trials, returns a N x 1 design
#'     matrix in which a set of weights (specified via the Mapping variable)
#'     are matched to the unique levels of X.
#'   }
#'
#' @return A design matrix
#' @examples
#' # Default is dummy coding
#' designCoding( 1:5 )
#' # Intercept, using mapping option
#' designCoding( 1:5, Levels = c(3,1,2,4,5), type = 'Intercept' )
#' # Example of effects coding with the levels reversed
#' # Note that here, mapping ranges from 0 - 4, not 1 - 5
#' designCoding( 1:5, Mapping = 4:0, type = 'Effects' )
#' # Coefficients
#' set.seed(500)
#' designCoding( 1:5, Mapping = rnorm(5), type = 'Coef' )
#'
#' @export
designCoding = function( X, Levels = NULL, Mapping=NULL,
                         type = 'Dummy' ) {

  # If necessary, extract the unique values of the variable
  if (length(Levels)==0) Levels = unique(X)

  # Determine the total number of levels in the variable
  K = length(Levels)

  # Dummy coding
  if (type=='Dummy') {

    # Create a design matrix
    out = matrix( 0, nrow = length(X), ncol = K-1 )

    # Create default mapping if necessary
    if (length(Mapping)==0) Mapping = 0:(K-1)

    for (k in 1:K) {

      if ( Mapping[k]!=0 ) out[X==Levels[k], Mapping[k]] = 1;

    }

  }

  # Effects coding
  if (type=='Effects') {

    # Create a design matrix
    out = matrix( 0, nrow = length(X), ncol = K-1 )

    # Create default mapping if necessary
    if (length(Mapping)==0) Mapping = 0:(K-1)

    for (k in 1:K) {

      if ( Mapping[k]==0 ) out[X==Levels[k]] = -1;
      if ( Mapping[k]!=0 ) out[X==Levels[k], Mapping[k]] = 1;

    }

  }

  # Intercept coding
  if (type=='Intercept') {

    # Create a design matrix
    out = matrix( 0, nrow = length(X), ncol = K )

    # Create default mapping if necessary
    if (length(Mapping)==0) Mapping = 1:K

    for (k in 1:K) {

      out[X==Levels[k], Mapping[k]] = 1;

    }

  }

  # Univariate Coefficient mapping
  if (type=='Coef') {

    if ( length(Mapping)==0) stop('Must provide weighting values')

    out = matrix( 0, nrow = length(X), ncol = 1 )

    for (k in 1:K) {
      out[X==Levels[k],1] = Mapping[k]
    }

  }

  out
}

# Lookup - 14
#' Softmax Function
#'
#' A generalization of the logistic function takes a K-dimensional vector of
#' arbitrary values and converts it to a K-dimensional vector of real values
#' in the range (0,1) that sum to 1. The function is also known as the
#' normalized exponential.
#'
#' The function can take either a vector of a matrix of values. If a matrix
#' the function is applied to each row of the matrix.
#'
#' @param x a vector of values from -Inf to Inf.
#'
#' @return a vector of values from 0 to 1 that sum to 1.
#' @examples
#' set.seed(3902)
#' ex = softmax( rnorm(5) )
#' sum( ex ) # Should equal 1
#' mat = matrix( rnorm(9), 3, 3 )
#' ex = softmax( mat )
#' rowSums( ex ) # Each row should sum to 1
#' @export

softmax = function(x) {

  # Vector case
  if ( is.vector( x ) ) {
    out = exp(x)/sum( exp(x) )
  }
  # Matrix case
  if ( is.matrix(x) ) {
    out = t( apply( x, 1, function(x) exp(x)/sum( exp(x) ) ) )
  }

  out
}

# Lookup - 15
#' Reverse of the Softmax Function
#'
#' A function that, given a vector of probabilities that sum to one, will
#' determine the closest values that could be passed to the softmax function
#' to produce those probabilities using R's optim function.
#'
#' @param y a vector of probabilities (positive, sum to 1).
#' @param init the starting values for the optim function. If empty, generates
#'   random values from a normal distribution.
#' @param restrict a logical vector indicating whether certain values of x
#'   should be fixed to 0.
#'
#' @return The output from the optim function.
#' @examples
#' set.seed(984)
#' input = rnorm(5)
#' sm = softmax( input )
#' output = reverseSoftmax( sm )
#' round(sm,3)
#' round(output$par,3)
#' @export

reverseSoftmax = function(y,init=NULL,restrict=NULL) {

  if (is.null(init)) init = rnorm(length(y));

  optim(init, function (x) { x[restrict] = 0; sum( ( y - softmax(x) )^2 )} )
}
