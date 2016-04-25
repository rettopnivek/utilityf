#----------------------------#
# Assorted utility functions #
#----------------------------#

# Package development
# library(devtools)
# library(roxygen2)

# Index
# Lookup - 01:  Standard Error of the Mean
# Lookup - 02:  Logit Function
# Lookup - 03:  Logistic Function
# Lookup - 04:  Softmax Function
# Lookup - 05:  Reverse of the Softmax Function
# Lookup - 06:  Error Function
# Lookup - 07:  Plotting Range Function
# Lookup - 08:  Generate a Blank Plot
# Lookup - 09:  Function to Draw Ellipses
# Lookup - 10:  Conversion between Degrees and Radians
# Lookup - 11:  Conversion to Cartesian Coordinates
# Lookup - 12:  Create a Design Matrix
# Lookup - 13:  Akaike's Information Criterion
# Lookup - 14:  Bayesion Information Criterion
# Lookup - 15:  Wrapper for Maximum Likelihood Estimation

### TO DO ###
# Move hdi function to different package
# Lookup - 07:  Highest Density Interval Estimator
# Change standard error extraction to S3 method
# Lookup - 10:  Calculate Standard Errors and Confidence Intervals

# Lookup - 01
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

# Lookup - 02
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

# Lookup - 03
#' Logistic Function
#'
#' This function applies the logistic function to a set of values
#'
#' @param x a set of values ( -Inf <= x <= Inf).
#' @return Returns a set of probabilities.
#' @export
#' @examples
#' logistic( c(0,-2.197225,2.1972255,-Inf,Inf) )

logistic = function(x) {
  1/(1+exp(-x))
}

# Lookup - 04
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

# Lookup - 05
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

# Lookup - 08
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

# Lookup - 09
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

# Lookup - 10
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

# Lookup - 11
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

# Lookup - 12
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

# Lookup - 15
#' Wrapper for Maximum Likelihood Estimation
#'
#' This function allows multiple runs of the optim function
#' with random starting values to better control for local
#' maxima/minima.
#'
#' @param dat the data to be fitted (formats can vary).
#' @param mle_fn a function to calculate the negative of the sum of the
#'   log-likelihood (must take three named inputs: 1. 'par' - a vector of
#'   parameters, 2. 'dat' - the data, 3. 'priors' - a variable for prior
#'   values (can be null)).
#' @param st_fn a function to generate a set of starting values
#'   (must take no parameters and output must match the input length
#'   for the mle_fn function).
#' @param grad_fn an optional function that returns a vector with the
#'   value of the first derivative for each parameter.
#' @param method a string giving the type of optimization routine that
#'   optim should use.
#' @param priors an optional variable to allow penalized maximum likelihood.
#' @param SE a logical value; if true, attempts to extract standard errors
#'   and confidence intervals for parameters.
#' @param emStop the number of attempts to find starting values that
#'   produce a defined likelihood function.
#' @param alpha the coverage interval for the confidence intervals around
#'   the parameters.
#' @param nRep the number of times to run the estimation routine.
#'   Running the routine multiple times with dispersed starting values
#'   increases the chances of finding the true maximum/minimum instead of
#'   local maxima/minima.
#' @param ... additional parameters for the list of control parameters.
#'   For instance, the upper limit for the number of iterations the
#'   optimization routine uses can be increased by 'maxit = 5000'.
#' @return Returns a list consisting of
#' \describe{
#'   \item{\code{param}}{a vector of the best-fitting parameter estimates.}
#'   \item{\code{logLik}}{the value for the sum of the log-likelihoods.}
#'   \item{\code{startVal}}{the vector of starting values used for the
#'     parameters.}
#'   \item{\code{convergenceCheck}}{a numerical code, when equal to 0
#'     indicates the model successfully converged.}
#'   \item{\code{hessianMatrix}}{the matrix for the hessian.}
#'   \item{\code{SE}}{a vector of standard errors for the parameters.}
#'   \item{\code{CI}}{a matrix with the lower and upper bounds for the
#'     confidence intervals around each parameter.}
#'   \item{\code{runTime}}{the time interval it took to estimate the
#'     parameters.}
#'   \item{\code{runTime}}{the original output given by \code{optim}.}
#' }
#' @export
#' @examples
#' mle_fn = function(par,dat,priors=NULL) sum(dnorm(dat,par[1],exp(par[2]),log=T)) # Likelihood function
#' st_fn = function() runif( 2, c(50,log(.2)),c(150,log(3)) ) # Starting values function
#' dat = rnorm( 100, 100, 15 ) # Generate data
#' results = MLE( dat, mle_fn, st_fn )
#' par = results$param; round( c(par[1],exp(par[2])) ) # The paramaters

MLE = function( dat, mle_fn, st_fn, grad_fn = NULL,
                method = 'Nelder-Mead', priors = NULL,
                SE = T, emStop = 20, alpha = .95, nRep = 5,
                ... ) {

  # Check that 'stats' package is installed
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("The 'stats' package is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  startTime = Sys.time()

  # Repeat data as a list nRep times
  datList = datList = rep( list(dat), nRep )

  # Run maximum likelihood wrapper nRep times
  results = lapply( datList, mleWrapper,
                    mle_fn = mle_fn,
                    st_fn = st_fn,
                    grad_fn = grad_fn,
                    priors = priors, SE = SE,
                    method = method,
                    emStop = emStop, ... )

  # Find results with the maximum log-likelihood value
  allLogLik = sapply( results, function(x) x[['logLik']] )

  # Determine for which iterations estimation failed
  estFail = which( is.na( allLogLik ) )
  if ( length( estFail ) == nRep )
    stop('Error: estimation failed for all iterations',call. = FALSE)

  # Find the iteration with the
  selBest = which.max( na.omit( allLogLik ) )
  # If there exist multiple iterations, take first one
  selBest = selBest[1]

  # Extract iteration with maximum log-likelihood
  resultsBest = results[[ selBest ]]

  # Initialize output
  output = list(
    param = NA,
    logLik = NA,
    startVal = NA,
    convergenceCheck = NA,
    hessianMatrix = NA,
    SE = NA,
    CI = NA,
    runTime = NA,
    optimOutput = NA
  )

  if ( is.list( resultsBest ) ) {

    output$param = resultsBest$par
    output$logLik = resultsBest$value
    output$startVal = resultsBest$startVal
    output$convergenceCheck = resultsBest$convergenceCheck

    if (SE) {

      output$hessianMatrix = resultsBest$hessianMatrix

      if ( sum( is.na( resultsBest$hessianMatrix ) ) == 0 ) {

        fisher_info = solve(-resultsBest$hessianMatrix)
        prop_sigma = sqrt(diag(fisher_info))

        ub = qnorm( (1 - alpha)/2 )
        lb = qnorm( (1 - alpha)/2, lower.tail = F )

        output$SE = prop_sigma
        output$CI = rbind( lb, ub )

      }

    }

    output$optimOutput = resultsBest$optimOutput

  }

  output$runTime = Sys.time() - startTime

  return( output )
}
