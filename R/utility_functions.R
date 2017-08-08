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
# Lookup - 13:  Extract Unique Levels from Combined Covariates
# Lookup - 14:  Create a Variable with Incremental Unit Intervals
# Lookup - 15:  Calculate Select Information Criterion Values
# Lookup - 16:  Draw a Violin Plot
# Lookup - 17:  Find the Mode
# Lookup - 18:  Draw Error Bars
# Lookup - 19:  Estimate Density for Individual Observations
# Lookup - 20:  Compute Category Proportions
# Lookup - 21:  Template for Documentation
# Lookup - 22:  Improved list
# Lookup - 23:  Calculate the First Four Moments
# Lookup - 24:  Generate Multi-line comments
# Lookup - 25:  Generate Custom Plot Axes

# Lookup - 01
#' Standard Error of the Mean
#'
#' This function calculates the standard error of the mean.
#'
#' @param x a vector of values.
#'
#' @return Returns the estimated standard error of the mean
#'   based on the values of x
#'
#' @examples
#' 1/sqrt(1000) # True value
#' sem( rnorm(1000) )
#'
#' @export

sem = function(x) {
  return( sd(x)/sqrt( length(x) ) )
}

# Lookup - 02
#' Logit Function
#'
#' This function calculates the logit (log of the odds) of a set of
#' probabilities.
#'
#' @param p a set of probabilities (0 <= p <= 1).
#'
#' @return Returns a set of values that now lie between -Inf and Inf.
#'
#' @examples
#' logit( c(.5,.1,.9,0,1) )
#'
#' @export

logit = function(p) {
  p[p<0 | p>1] = NA
  return( log(p/(1-p)) )
}

# Lookup - 03
#' Logistic Function
#'
#' This function applies the logistic function to a set of values
#'
#' @param x a set of values ( -Inf <= x <= Inf).
#'
#' @return Returns a set of probabilities.
#'
#' @examples
#' logistic( c(0,-2.197225,2.1972255,-Inf,Inf) )
#'
#' @export

logistic = function(x) {
  return( 1/(1+exp(-x)) )
}

# Lookup - 04
#' Softmax Function
#'
#' A generalization of the logistic function takes a K-dimensional
#' vector of arbitrary values and converts it to a K-dimensional
#' vector of real values in the range (0,1) that sum to 1. The
#' function is also known as the normalized exponential.
#'
#' The function can take either a vector or a matrix of values.
#' If a matrix is passed in, the function is applied to each row
#' of the matrix.
#'
#' @param x a vector of values from -Inf to Inf.
#'
#' @return a vector of values from 0 to 1 that sum to 1.
#'
#' @examples
#' set.seed(3902)
#' ex = softmax( rnorm(5) )
#' sum( ex ) # Should equal 1
#' mat = matrix( rnorm(9), 3, 3 )
#' ex = softmax( mat )
#' rowSums( ex ) # Each row should sum to 1
#'
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

  return( out )
}

# Lookup - 05
#' Reverse of the Softmax Function
#'
#' A function that, given a vector of probabilities that sum to one,
#' will determine the closest values that could be passed to the
#' softmax function to produce those probabilities using the
#' \code{\link[stats]{optim}} function.
#'
#' @param y a vector of probabilities (positive, sum to 1).
#' @param init the starting values for the optim function. If empty,
#'   generates random values from a normal distribution.
#' @param restrict a logical vector indicating whether certain values
#'   of x should be fixed to 0.
#'
#' @return The output from the optim function.
#'
#' @examples
#' set.seed(984)
#' input = rnorm(5)
#' sm = softmax( input )
#' output = reverseSoftmax( sm )
#' round(sm,3)
#' round(output$par,3)
#'
#' @export

reverseSoftmax = function(y,init=NULL,restrict=NULL) {

  if (is.null(init)) init = rnorm(length(y));

  # Least-squares difference
  f = function (x) { x[restrict] = 0; sum( ( y - softmax(x) )^2 )}

  # Minimize function
  out = optim(init, f )

  return( out )
}

# Lookup - 06
#' Error Function
#'
#' Calculates the error function.
#'
#' @param x a vector of values on the real number line.
#'
#' @return a vector of transformed values based on the error function.
#'
#' @examples
#' plot( c(-1,1), c(-3,3), type='n', xlab='x', ylab='erf(x)' )
#' x = seq(-3,3,length=100)
#' lines(x,erf(x))
#'
#' @export

erf = function(x) {
  return( 2*pnorm( x*sqrt(2), 0, 1 ) - 1 )
}

# Lookup - 07
#' Plotting Range Function
#'
#' This function determines the lower and upper boundaries for
#' a vector of values given a specified sub-interval for
#' plotting purposes.
#'
#' @param int the desired sub-interval over the range.
#' @param dat the vector of values over which to determine the
#'   boundaries.
#'
#' @return Returns a lower and upper boundary, multiples of the
#'   specified sub-interval.
#'
#' @examples
#' d = density( rnorm(100) )
#' x.limit = lowerUpper( .5, d$x )
#' y.limit = lowerUpper( .1, d$y )
#' plot( x.limit, y.limit, type='n', xlab='Values', ylab='Density' )
#' lines(d$x,d$y)
#'
#' @export

lowerUpper = function(int,dat) {

  ll = int*floor(min(dat)/int)
  ll = c(ll,int*ceiling(max(dat)/int))

  return( ll )
}

# Lookup - 08
#' Generate a Blank Plot
#'
#' This function generates a completely blank plot.
#'
#' @param xDim the lower and upper boundaries for the
#'   x-axis.
#' @param yDim the lower and upper boundaries for the
#'   y-axis.
#'
#' @return A empty plot.
#'
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
#'   for the ellipse.
#' @param Xc the center on the x-axis for the ellipse.
#' @param Yc the center on the y-axis for the ellipse.
#' @param deg if TRUE, assumes k is in degrees; otherwise,
#'   k is assumed to be in radians.
#' @param lngth the number of points to use when plotting
#'   the ellipse.
#' @param draw if TRUE, draws a polygon in the shape of the
#'   ellipse on an existing plot.
#' @param ret if TRUE, returns a matrix with the x and y
#'   values for the ellipse.
#' @param ... additional parameters for the polygon function.
#'
#' @return If \code{ret} is set to TRUE, returns the x and y-axis
#'   values that were passed into the \code{\link[graphics]{polygon}}
#'   function.
#'
#' @examples
#' plot( c(-1,1), c(-1,1), type='n', xlab='X-axis', ylab='Y-axis' )
#' drawEllipse( .2, .2 )
#' ex = drawEllipse( .5, .2, Xc = -.5, Yc = .5, k = 20, ret = T, col = 'green')
#' drawEllipse( .05, .3, Xc = .5, Yc = -.5, k = 290, col = 'red')
#'
#' @export

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
#'
#' @examples
#' r = degreesRadians( degrees = 45 )
#' print( r )
#' d = degreesRadians( radians = .7853982 )
#' print( d )
#'
#' @export

degreesRadians = function( degrees = NULL, radians = NULL ) {

  out = NULL

  if ( length( degrees )>0) {
    out = degrees*(pi/180)
  }

  if ( length( radians )>0) {
    out = radians/(pi/180)
  }

  return( out )
}

# Lookup - 11
#' Conversion to Cartesian Coordinates
#'
#' A function to use the magnitude and angle of a 2D vector
#' and convert them to cartesian coordinates (assumes the
#' vector starts at [0,0]).
#'
#' @param H a magnitude.
#' @param A an angle.
#' @param degrees if true, the angle is assumed to be in
#'   degrees.
#'
#' @return A set of cartesian coordinates.
#'
#' @examples
#' cart = convertMagnitudeAngle( 2, 45 )
#' print( cart )
#'
#' @export

convertMagnitudeAngle = function( H, A, degrees = T) {

  # Convert to radians
  if (degrees) A = degreesRadians( degrees = A )

  a = cos(A)*H; b = sin(A)*H

  return( c(a,b) )
}

# Lookup - 12
#' Create Design Matrix
#'
#' A function that translates a vector of categories into different
#' types of design matrices.
#'
#' @param X a vector of categories
#' @param Levels the unique values of X. If left unspecified, the
#'   function attemps to extract, but the mappings will be based on
#'   order of appearance.
#' @param Mapping For \code{Dummy} and \code{Effects} options,
#'   provides an additional way to indicates which unique values of
#'   X should be mapped to which columns of the design matrix.
#'   For the \code{Coef} option, provides the corresponding weight
#'   values for the unique values of X.
#' @param type Current options include...
#'   \itemize{
#'     \item \code{Dummy} -> Given K levels and N trials, returns a
#'     N x (K-1) design matrix using dummy coding (i.e. 0 and 1,
#'     with one variable). Useful for simple effects.
#'     \item \code{Effects} -> Given K levels and N trials, returns a
#'     N x (K-1) design matrix using effects coding (i.e. -1 and 1,
#'     with one variable denoted solely by -1). Useful for comparisons
#'     against the grand mean.
#'     \item \code{Intercept} -> Given K levels and N trials, returns a N x K
#'     design matrix in which each unique level has its own column
#'     (i.e. K unique intercepts).
#'     \item \code{Coef} -> Given K levels and N trials, returns a N x 1
#'     design matrix in which a set of weights (specified via the
#'     \code{Mapping} variable) are matched to the unique levels of X.
#'   }
#'
#' @return A design matrix.
#'
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

  return( out )
}

# Lookup - 13
#' Extract Unique Levels from Combined Covariates
#'
#' This function creates a single variable with a set of unique
#' levels based on a set of covariates, the number of all possible
#' combinations for each of the covariates
#'
#' @param mat a matrix of the covariates of interest.
#'
#' @return Returns a vector indicating the corresponding combination of
#'   levels for each observation.
#'
#' @examples
#' # Create an example matrix of covariates with differing levels
#' mat = cbind( rep( c(0,1), 3 ), rep( c(0,1,0), each = 2 ) )
#' ex = covCreate( mat ); print( ex )
#'
#' @export

covCreate = function(mat) {

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

  return( out )
}

# Lookup - 14
#' Create a Variable with Incremental Unit Intervals
#'
#' A function to take a discrete variable and generate a new variable
#' with the same number of levels whose values increase in unit
#' increments. This is useful, for instance, if you have a variable
#' encoding irregular subject ID numbers that you would like to use
#' for indexing purposes.
#'
#' @param x a vector of values.
#'
#' @return Returns a vector of matching length to \code{x} with
#'   incremental unit intervals corresponding to each level of the
#'   original variable.
#'
#' @examples
#' x = c( 2320, 1038, 3010, 7503 )
#' print( createIncrement(x) )
#'
#' @export

createIncrement = function(x) {

  # Determine the total number of original values
  curVal = sort( unique( x ) )
  newVal = 1:length(curVal) # Create regular increments
  newX = numeric( length( x ) )
  for (i in 1:length(curVal)) {
    sel = x == curVal[i]
    newX[sel] = newVal[i]
  }

  return( newX )
}


# Lookup - 15
#' Calculate Select Information Criterion Values
#'
#' A function that calculates either Akaike's Information Criterion (AIC)
#' with a correction or the Bayesian Information Criterion (BIC).
#'
#' @param logLik a log-likelihood value.
#' @param k the number of free parameters for the model.
#' @param n the number of observations in the sample.
#' @param type indicates whether the 'AICc' or 'BIC' value should be returned.
#'
#' @details Given a summed log-likelihood \eqn{L} from a model and \eqn{K} free
#'   parameters, the AIC is \deqn{ 2K - 2L. } A correction for finite samples is
#'   recommended (e.g. Burnham & Anderson, 2002), and for \eqn{N} observations
#'   the new equation is \deqn{ 2K - 2L + \frac{2K(K+1)}{N+K+1}. } The formula
#'   for the BIC is \deqn{ log(N)K - 2L. } For both criterions, models with
#'   smaller values are to be preferred.
#'
#' @references
#' Akaike, H. (1973). Information theory and an extension of the maximum likelihood
#'   principle. In B. N. Petrov & F. Caski (Eds.), Proceedings of the Second
#'   International Symposium on Information Theory (pp. 267-281). Budapest:Akademiai
#'   Kiado.
#'
#' Burnham, K. P., & Anderson, D. R. (2002). Model selection and multimodel inference:
#'   A practical information-theoretic approach. New York: Springer-Verlag.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. Annals of Statistics, 6,
#'   461-464.
#'
#' @return A value for either the AICc or the BIC.
#'
#' @examples
#' N = 100; K = 2
#' x = rnorm( 100, 1, 2 )
#' m1 = sum( dnorm( x, 0, 1, log = T ) )
#' m2 = sum( dnorm( x, 1, 2, log = T ) )
#' # AIC values comparing the two models
#' print( round( infoCrit( c(m1,m2), K, N ), 2 ) )
#' # BIC values comparing the two models
#' print( round( infoCrit( c(m1,m2), K, N ), 2 ), type = 'BIC' )
#'
#' @export

infoCrit = function( logLik, k, n, type = 'AICc' ) {

  # Initialize output
  out = NA

  if ( type == 'AICc' ) out = AICc_val( logLik, k, n )
  if ( type == 'BIC' ) out = BIC_val( logLik, k, n )

  return( out )
}

# Lookup - 16
#' Draw a Violin Plot
#'
#' Adds a violin plot to an already existing figure.
#'
#' @param x a vector of values used to create the
#'   violin plot, or a named list where \code{X} is
#'   the x-axis values and \code{y} is the associated
#'   density.
#' @param pos the x-axis position at which to draw the
#'   violin plot
#' @param scaleH the maximum half-width of the violin plot.
#' @param est a logical value; if true, the density is
#'   empirically estimated, otherwise the density is
#'   extracted from the list \code{x}.
#' @param interval allows for extra options governing ways
#'   to restrict the violin plot to a  desired sub-interval.
#'   \describe{
#'     \item{\code{crit}}{a critical value or pair of values
#'       that determines the sub-interval to plot over.}
#'     \item{\code{type}}{when \code{crit} takes on a single
#'       value, 'greater' restricts the sub-interval to values
#'       larger than \code{crit},'less' to values below
#'       \code{crit}.}
#'     \item{\code{out}}{a logical value; if true, the
#'       function returns the proportion of values above,
#'       below, or within \code{crit}.}
#'   }
#' @param ... additional plotting parameters for
#'   \code{\link[graphics]{polygon}}.
#'
#' @return If \code{out} is set to TRUE, returns the proportion of
#'   values in \code{x} above, below, or within the interval defined
#'   by \code{crit}.
#'
#' @examples
#' x = rnorm(100)
#' blankPlot(yDim=c(-6,6))
#' crit = mean(x)
#' violinPlot( x, .5, interval = list( crit = crit ), col = 'grey',
#'   border = NA )
#' violinPlot( x, .5, lwd = 2 )
#'
#' @export

violinPlot = function( x, pos, scaleH = .5, est = T,
                       interval = list( ), ... ) {

  # Set default options for graphing a specific interval
  if ( length( interval$crit ) == 0 ) crit = NULL else
    crit = interval$crit
  if ( length( interval$type ) == 0 ) type = 'greater' else
    type = interval$type
  if ( length( interval$out ) == 0 ) out = F else
    out = interval$out

  # Calculate density
  if (est) den = density( x ) else { den = x; x = den$x }

  # Restrict to a specific interval
  if ( length(crit) == 2 ) {
    sel = den$x >= crit[1] & den$x <= crit[2]
    PostProb = sum( x >= crit[1] & x <= crit[2] )/length(x)
  }
  # Restrict to a specific half
  if ( length(crit) == 1 ) {
    if (type=='greater') {
      sel = den$x > crit
      PostProb = sum( x > crit )/length(x)
    }
    if (type=='less') {
      sel = den$x < crit
      PostProb = sum( x < crit )/length(x)
    }
  }
  # No interval
  if ( length(crit) == 0 ) {
    sel = rep( T, length( den$x ) )
    PostProb = 1
  }

  den$y = den$y/max(den$y); den$y = den$y*scaleH;
  xa = c( -den$y[sel], rev(den$y[sel]) ) + pos
  ya = c( den$x[sel], rev(den$x[sel]) )
  polygon( xa, ya, ... )

  if ( out ) return( PostProb )
}

# Lookup - 17
#' Find the Mode
#'
#' This function determines the mode for a vector.
#'
#' @param x a vector (can be numeric or a set of strings).
#' @param discrete a logical value; if true, an empirical
#'   estimate of the probability mass function is used to
#'   determine the mode. If false, an empirical estimate of
#'   the probability density function is used instead.
#' @param ... additional options for the \code{density}
#'   function.
#'
#' @return Returns the estimated mode for the vector.
#'
#' @examples
#' x = rnorm(100)
#' findMode(x)
#' x = rbinom(10,20,.5)
#' findMode(x,discrete=T)
#' x = c( 'cat','dog','cat','mouse')
#' findMode(x,discrete=T)
#'
#' @export


findMode = function( x, discrete = F, ... ) {

  if ( discrete ) {

    epmf = table( x )/length( x )
    md = names( which( epmf == max( epmf ) ) )
    if ( mode( x ) == 'numeric' ) md = as.numeric( md )

  } else {
    epdf = density(x,...);
    md = epdf$x[ which( max(epdf$y) == epdf$y ) ]
  }

  return( md )
}

# Lookup - 18
#' Draw Error Bars
#'
#' Adds error bars to an already existing plot.
#'
#' @param pos a single value or a vector of N values, indicating
#'   the position(s) at which an error bar should be drawn.
#' @param limits a vector of 2 values or a 2 x N matrix giving the
#'   lower and upper limits, respectively, of the uncertainty
#'   intervals.
#' @param flip a logical value. If true, bars are drawn horizontally
#'   instead of vertically. In this case \code{pos} denotes the
#'   position(s) on the y-axis. Otherwise, \code{pos} denotes the
#'   position(s) on the x-axis.
#' @param ... additional plotting parameters for the
#'   \code{\link[graphics]{arrows}} function.
#'
#' @examples
#' # Example of 95% versus 68% intervals for standard normal
#' plot( c(0,1),c(-6,6),type='n',xaxt='n',xlab=' ',ylab='z-scores' )
#' pos = c(.25,.75)
#' limits = cbind( qnorm( c(.025,.975) ), qnorm( c(.16,.84) ) )
#' errorBars( pos, limits, length = .05, lwd = 2 )
#' text( pos, limits[2,], c("95%","68%"), pos = 3 )
#' abline(h=0)
#'
#' @export

errorBars = function( pos, limits, flip = F, ... ) {

  if ( flip ) {
    if ( is.matrix( limits ) ) {
      arrows( limits[1,], pos,
              limits[2,], pos, code = 3, angle = 90, ... )
    } else {
      arrows( limits[1], pos,
              limits[2], pos, code = 3, angle = 90, ... )
    }
  } else {
    if ( is.matrix( limits ) ) {
      arrows( pos, limits[1,],
              pos, limits[2,], code = 3, angle = 90, ... )
    } else {
      arrows( pos, limits[1],
              pos, limits[2], code = 3, angle = 90, ... )
    }
  }

}

# Lookup - 19
#' Estimate Density for Individual Observations
#'
#' Given a vector of values, computes the empirical density
#' for each observation.
#'
#' @param x a numeric vector.
#' @param ... additional parameters for the
#'   \code{\link[stats]{density}} function.
#'
#' @return A list with...
#'   \itemize{
#'     \item x - the sorted values for the original input.
#'     \item y - the associated empirical densities.
#'   }
#'
#' @examples
#' plot(c(-4,4),c(0,.5),type='n',ylab='Density',xlab='z-scores')
#' x = rnorm( 100 )
#' dp = densityPoints( x )
#' points( dp$x, dp$y, pch = 19 )
#'
#' @export

densityPoints = function( x, ... ) {

  ed = density( x, ... )
  af = approxfun( ed )
  y = af( sort( x ) )

  return( list( x = sort( x ), y = y ) )
}

# Lookup - 20
#' Compute Category Proportions
#'
#' Given a set of posssible values, computes the proportion of
#' values that actually occurred.
#'
#' @param x a vector of values.
#' @param val a character vector with the set of possible values
#'   that \code{x} can take on.
#'
#' @return Returns a vector of proportions that sum to 1.
#'
#' @examples
#' x = rbinom( 100, 1, .7 )
#' catProp( x, c("0","1") )
#'
#' # Multiple categories, some categories don't occur
#' u = runif(100)
#' x = rep( 'dog', 100 )
#' x[ u < .33 ] = 'cat'
#' catProp( x, c("cat","dog","bird" ) )
#'
#' @export

catProp = function( x, val ) {

  # Initialize output
  out = numeric( length( val ) )
  names( out ) = val

  # Tally frequencies for each choice
  freq = table( x )
  # Compute associated proportions
  prp = freq/sum( freq )

  # Extract values
  out[ names( prp ) ] = as.numeric( prp )

  return( out )
}

# Lookup - 21
#' Template for Documentation
#'
#' Prints a template for an informal way to document R functions
#' to the console.
#'
#' @return Prints the template to the console window, where it
#'   can be copied and modified by the user.
#'
#' @examples
#' quickDocTemplate()
#'
#' @export

quickDocTemplate = function() {
  cat( "  # Purpose:", "\n" )
  cat( "  # ...", "\n" )
  cat( "  # Arguments:", "\n" )
  cat( "  # ...", "\n" )
  cat( "  # Returns:", "\n" )
  cat( "  # ...", "\n" )
}

# Lookup - 22

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

#' Improved List of Objects
#'
#' A function created by Dirk Eddelbuettel that lists the objects
#' currently in the workspace along with pertinent information.
#'
#' @return A list of objects in the workspace, along with their
#' memory usage and dimensions.
#'
#' @section References:
#' Eddelbuettel, D. (2009). Retrieved from http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
#'
#' @examples
#' x = rnorm( 100 )
#' lsos()
#'
#' @export

lsos = function( ..., n = 10 ) {
  .ls.objects( ..., order.by = "Size", decreasing = TRUE,
               head = TRUE, n = n )
}

# Lookup - 23
#' Calculate the First Four Moments
#'
#' A function to compute the first four moments (i.e., the mean,
#' variance, skew, and kurtosis) for a vector of data.
#'
#' @param x a numerical vector of values.
#'
#' @return A list consisting of...
#'   \itemize{
#'     \item \code{n}: The number of observations.
#'     \item \code{mean}: The mean of the data (The 1st moment).
#'     \item \code{var}: The variance of the data (The 2nd moment).
#'     \item \code{sd}: The standard deviation of the data.
#'     \item \code{skew}: The skewness of the data (The 3rd moment).
#'     \item \code{kurtosis}: The kurtosis of the data (The
#'       fourth moment).
#'   }
#'
#' @section References:
#' Terriberry, T. B. (2007). Computing Higher-Order Moments Online.
#    Retrieved from https://people.xiph.org/~tterribe/notes/homs.html
#'
#' @examples
#' x = rnorm( 100 )
#' higherOrderMoments( x )
#'
#' @export

higherOrderMoments = function( x ) {

  n = 0 # Sample size
  M = 0 # Mean
  M2 = 0 # 2nd moment
  M3 = 0 # 3rd moment
  M4 = 0 # 4th moment

  # Function to raise a value to a power
  pow = function( v, a ) return( v^a )

  # Loop over observations
  for ( k in 1:length( x ) ) {

    n1 = n; n = n + 1;
    term1 = x[k] - M;
    term2 = term1 / n;
    term3 = term1 * term2 * n1;
    # Update mean
    M = M + term2;
    # Update fourth moment
    M4 = M4 + term3 * pow( term2, 2.0 ) * ( n*n - 3.0 * n + 3.0 ) +
      6.0 * pow( term2, 2.0 ) * M2 - 4.0 * term2 * M3;
    # Update third moment
    M3 = M3 + term3 * term2 * (n - 2.0) - 3.0 * term2 * M2;
    # Update second moment
    M2 = M2 + term3;
  }

  # Compute descriptive statistics from moments
  out = c(
    n = n,
    mean = M,
    var = M2 / ( n - 1.0 ),
    sd = sqrt( M2 / (n - 1.0 ) ),
    skewness = sqrt( n ) * M3/ pow(M2, 1.5),
    kurtosis = n*M4 / (M2*M2) - 3.0
  )

  return( out )
}

# Lookup - 24
#' Generate Multi-line comments
#'
#' A function that allows the user to write multi-line
#' comments that are not run.
#'
#' @param string A character string with the desired comments.
#'
#' @examples
#' b("
#' First line of comments
#' Second line of comments
#' ")
#'
#' @export

b = function( string ) {
  # Silent return
  if (0) {
    string
  }
}

# Lookup - 25
#' Generate Custom Plot Axes
#'
#' A function that, based on the x and y-axis boundaries,
#' draws lines and labels for a specified set of axes.
#'
#' @param xl a vector with the lower and upper boundaries for the x-axis.
#' @param yl a vector with the lower and upper boundaries for the y-axis.
#' @param pos a vector indicating at which sides an axis should be drawn.
#'   The possible sides are:
#'   \enumerate{
#'     \item Bottom (codes = \code{'Bottom'}, \code{'B'}, \code{'b'},
#'       \code{'1'}, \code{1});
#'     \item Left (codes = \code{'Left'}, \code{'L'}, \code{'l'},
#'       \code{'2'}, \code{2});
#'     \item Top (codes = \code{'Top'}, \code{'T'}, \code{'t'},
#'       \code{'3'}, \code{3});
#'     \item Right (codes = \code{'Right'}, \code{'R'}, \code{'r'},
#'       \code{'4'}, \code{4});
#'   }
#'   The number of elements in the vector tell the function how many axes
#'   to draw.
#' @param lnSz the width of the axis lines.
#' @param type the type of axis line to draw. \code{1} or \code{'extend'}
#'   generates axis lines that extend past the boundaries, while \code{2}
#'   or \code{'truncated'} truncates the axis lines at the boundaries.
#' @param label the labels for the axes. Must match the length of variable
#'   \code{pos} else a warning is generated and no labels are added.
#' @param lbPos the line positions of the labels for the axes.
#' @param lbSz the font-size of the labels for the axes.
#' @param inc the increments for the axis values. Must match the length of
#'   the variable \code{pos} else a warning is generated and no axis values
#'   are added. Set to \code{0} to suppress the adding of axis values.
#' @param axisLabel An optional list matching in length to the variable
#'   \code{pos}, with each element consisting either of a list of the tick
#'   positions and their labels for a given axis, or \code{NULL} to suppress
#'   adding values to that axis.
#' @param axSz the size of the font for the axis values.
#' @param axPos the line position of the axes.
#' @param prec the number of decimal values to report for the axis values.
#'
#' @examples
#' layout( cbind( 1, 2 ) ) # 2 plotting panes
#' xl = c(0,1); yl = c(0,1); blankPlot(xl,yl) # Empty plot
#' customAxes( xl, yl, label = c( 'X', 'Y' ), lnSz = 2,
#'   inc = c(.25,.25), axSz = .75 )
#' xl = c(.5,2.5); yl = c(-1,1); blankPlot(xl,yl); abline( h = 0 )
#' customAxes( xl, yl, pos = c(1,4), label = c( 'Conditions', 'Values' ),
#'   lbSz = 1, lbPos = c(1,1),
#'   axisLabel = list( list( c(1,2), c('A','B') ), NULL ),
#'   inc = c(0,.5), axSz = .75 )
#'
#' @export
customAxes = function( xl, yl,
                       pos = c( 1, 2 ),
                       lnSz = 2,
                       type = 2,
                       label = NULL,
                       lbPos = c( 2.5, 2.5, 1, 2.5 ),
                       lbSz = 1.5,
                       inc = c(0,0,0,0),
                       axisLabel = NULL,
                       axSz = 1.25,
                       axPos = c( -1, -1, -1, -1 ),
                       prec = 2 ) {

  # Check whether labels for axes should be added
  labCheck = F
  if ( !is.null( label ) ) {
    if ( length( label ) != length( pos ) ) {
      warning( "Number of labels does not match number of specified axes" )
    } else labCheck = T
  }

  alCheck = F
  if ( !is.null( axisLabel ) ) {
    if ( !is.list( axisLabel ) ) {
      warning( "A list of axis labels is needed" )
    } else if ( length( axisLabel ) != length( pos ) ) {
      warning( "Number of axis labels does not match number of specified axes" )
    } else alCheck = T
  }

  for ( i in 1:length( pos ) ) {

    # Bottom axis
    if ( pos[i] == 1 | pos[i] == '1' | pos[i] == 'B' |
         pos[i] == 'b' | pos[i] == 'bottom' |
         pos[i] == 'Bottom' ) {

      # Draw axis line
      if ( type == 1 | type == 'extend' )
        abline( h = yl[1], lwd = lnSz )
      if ( type == 2 | type == 'truncate' )
        segments( xl[1], yl[1], xl[2], yl[1], lwd = lnSz )

      # Add label to axis
      if ( labCheck ) {
        mtext( label[i], side = 1, cex = lbSz, line = lbPos[i] )
      }

      # Draw axis tick labels

      # Increments
      if ( inc[i] > 0 ) {
        ai = seq( xl[1], xl[2], inc[i] )
        ai = round( ai, prec )
        axis( 1, ai, tick = F, line = axPos[i],
              cex.axis = axSz )
      }

      # String labels
      if ( alCheck ) {
        if ( !is.null( axisLabel[[i]] ) ) {
          axis( 1, axisLabel[[i]][[1]], axisLabel[[i]][[2]],
                tick = F, line = axPos[i],
                cex.axis = axSz )
        }
      }
    }

    # Left axis
    if ( pos[i] == 2 | pos[i] == '2' | pos[i] == 'L' |
         pos[i] == 'l' | pos[i] == 'left' |
         pos[i] == 'Left' ) {

      # Draw axis line
      if ( type == 1 | type == 'extend' )
        abline( v = xl[1], lwd = lnSz )
      if ( type == 2 | type == 'truncate' )
        segments( xl[1], yl[1], xl[1], yl[2], lwd = lnSz )

      # Add axis label
      if ( labCheck ) {
        mtext( label[i], side = 2, cex = lbSz, line = lbPos[i] )
      }

      # Draw axis tick labels

      # Increments
      if ( inc[i] > 0 ) {
        ai = seq( yl[1], yl[2], inc[i] )
        ai = round( ai, prec )
        axis( 2, ai, tick = F, line = axPos[i],
              cex.axis = axSz )
      }

      # String labels
      if ( alCheck ) {
        if ( !is.null( axisLabel[[i]] ) ) {
          axis( 2, axisLabel[[i]][[1]], axisLabel[[i]][[2]],
                tick = F, line = axPos[i],
                cex.axis = axSz )
        }
      }
    }

    # Top axis
    if ( pos[i] == 3 | pos[i] == '3' | pos[i] == 'T' |
         pos[i] == 't' | pos[i] == 'top' |
         pos[i] == 'Top' ) {

      # Draw axis line
      if ( type == 1 | type == 'extend' )
        abline( h = yl[2], lwd = lnSz )
      if ( type == 2 | type == 'truncate' )
        segments( xl[1], yl[2], xl[2], yl[2], lwd = lnSz )

      # Add axis label
      if ( labCheck ) {
        mtext( label[i], side = 3, cex = lbSz, line = lbPos[i] )
      }

      # Draw axis tick labels

      # Increments
      if ( inc[i] > 0 ) {
        ai = seq( xl[1], xl[2], inc[i] )
        ai = round( ai, prec )
        axis( 3, ai, tick = F, line = axPos[i],
              cex.axis = axSz )
      }

      # String labels
      if ( alCheck ) {
        if ( !is.null( axisLabel[[i]] ) ) {
          axis( 3, axisLabel[[i]][[1]], axisLabel[[i]][[2]],
                tick = F, line = axPos[i],
                cex.axis = axSz )
        }
      }


    }

    # Right axis
    if ( pos[i] == 4 | pos[i] == '4' | pos[i] == 'R' |
         pos[i] == 'r' | pos[i] == 'right' |
         pos[i] == 'Right' ) {

      # Draw axis line
      if ( type == 1 | type == 'extend' )
        abline( v = xl[2], lwd = lnSz )
      if ( type == 2 | type == 'truncate' )
        segments( xl[2], yl[1], xl[2], yl[2], lwd = lnSz )

      # Add axis label
      if ( labCheck ) {
        mtext( label[i], side = 4, cex = lbSz, line = lbPos[i] )
      }

      # Draw axis tick labels

      # Increments
      if ( inc[i] > 0 ) {
        ai = seq( yl[1], yl[2], inc[i] )
        ai = round( ai, prec )
        axis( 4, ai, tick = F, line = axPos[i],
              cex.axis = axSz )
      }

      # String labels
      if ( alCheck ) {
        if ( !is.null( axisLabel[[i]] ) ) {
          axis( 1, axisLabel[[i]][[1]], axisLabel[[i]][[2]],
                tick = F, line = axPos[i],
                cex.axis = axSz )
        }
      }

    }

  }

}
