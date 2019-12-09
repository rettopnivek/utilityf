#----------------------------#
# Assorted utility functions #
#----------------------------#

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
# Lookup - 22:  Improved List of Objects
# Lookup - 23:  Calculate the First Four Moments
# Lookup - 24:  Generate Multi-line comments
# Lookup - 25:  Generate Custom Plot Axes
# Lookup - 26:  Relative Weights for Information Criterion Values
# Lookup - 27:  Simple Bootstrap Function
# Lookup - 28:  Template for Project Headers
# Lookup - 29:  Add Horizontal Lines to a Plot
# Lookup - 30:  Add Vertical Lines to a Plot
# Lookup - 31:  Fit a Basic Piecewise Regresssion
# Lookup - 32:  Add Basic Color Map to a Plot
# Lookup - 33:  Raise a Value to a Power
# Lookup - 34:  Custom Function to Standardize a Variable
# Lookup - 35:  Functions to Evaluate Run Times
# Lookup - 36:  Function to Identify NA Values by Row
# Lookup - 37:  Transform Hit/False Alarm Rates into SDT Parameters
# Lookup - 38:  Estimate Parameters for a Two-Boundary Unbiased Wiener Process
# Lookup - 39:  Plot a Heatmap for a Correlation Matrix
# Lookup - 40:  Template for Standard Plotting Code
# Lookup - 41:  Print a Nicely Formatted Table
# Lookup - 42:  Estimate or Format p-values from Monte Carlo Samples
# Lookup - 43:  Colorblind Friendly Palette

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
#' Coding Schemes for Nominal Variables
#'
#' Creates a numeric coding scheme for a nominal variable
#' with two or more categories.
#'
#' @param x a vector of levels.
#' @param type the type of coding to use, either 'Dummy',
#' 'Effects', or 'Intercept'.
#' @param levels an optional value used to specify the
#'   reference group with dummy, simple, and effects
#'   coding schemes, or more generally, a vector with
#'   the unique levels used to map a desired order when
#'   creating the design matrix.
#' @param label an optional character string giving the
#'   label for the nominal variable.
#' @param weights an optional vector of weights to be
#'   assigned to each level of the nominal variable
#'   when applying the coefficient coding scheme.
#'
#' @details
#'
#' With dummy coding, each additional category in a nominal
#' variable is compared against a reference category. For
#' K categories, K - 1 dummy variables are created, coded as 1
#' for the presence of a category and 0 otherwise. The intercept
#' term is interpreted as the cell mean for the reference category.
#'
#' With simple coding, each additional category in a nominal
#' variable is compared against a reference category. For
#' K categories, K - 1 dummy variables are created, coded as
#' (K-1)/K for the presence of a category, -1/K otherwise.
#' The intercept term is interpreted as the grand mean
#' (the mean of the cell means).
#'
#' With effects coding (also known as deviation coding),
#' each additional category in a nominal variable is
#' compared against the grand mean. For K categories, K - 1
#' dummy variables are created, coded as 1 for the presence
#' of a category, -1 for the presence of the reference
#' category, and 0 otherwise. The intercept term is
#' interpreted as the grand mean.
#'
#' With intercept coding, a separate dummy variable is specified
#' for each category in the nominal variable, coded as 1 when
#' the category is present and 0 otherwise. This coding scheme
#' requires that the model have no intercept term. Instead,
#' predictions for the dependent variable for each category are
#' estimated separately. Therefore, for K categories, K dummy
#' variables are created.
#'
#' With coefficient coding, a single dummy variable is specified,
#' and a desired weight is assigned for each level of the nominal
#' variable. This is useful for testing hypothesized ordered
#' relationships.
#'
#' @return Given K levels for the inputted variable, returns a matrix
#' with K - 1 columns (dummy, simple, and effects coding schemes),
#' K columns (intercept coding schemes), or 1 column (coefficient
#' coding schemes) with the new numerical coding scheme.
#'
#' @examples
#' x = rep( c('low','med','high'), each = 2 ) # 3 levels
#' data.frame( x, nominalCoding( x, type = 'Dummy' ) )
#' data.frame( x, nominalCoding( x, type = 'Simple' ) )
#' data.frame( x, nominalCoding( x, type = 'Effects' ) )
#' data.frame( x, nominalCoding( x, type = 'Intercept' ) )
#' data.frame( x, nominalCoding( x, type = 'Coefficient', weights = rnorm(3) ) )
#'
#' @export
nominalCoding = function( x, type = 'Effects',
                          levels = NULL, label = NULL,
                          weights = NULL ) {

  # Extract variable name
  if ( is.null( label ) ) {
    nm = deparse( substitute(x) )
  } else {
    nm = label
  }

  # Extract dimensions
  n = length( x ) # Number of values
  k = length( unique( x ) ) # Number of levels

  # Extract unique levels
  if ( is.null( levels ) ) {
    levels = unique( x )
  } else {

    # Check that values in variable match inputted levels
    if ( !all( levels %in% x ) ) {
      stop( 'Levels must match values in variable',
            call. = FALSE )
    }

    # If a single level is provided as a reference
    if ( length( levels ) == 1 ) {

      # Reorder extracted levels so that reference
      # group is first
      tmp = unique( x )
      tmp = tmp[ tmp != levels ]
      levels = c( levels, tmp )

    } else {

      # If fewer categories than actual number of levels
      # are provided
      if ( length( levels ) < k ) {

        tmp = unique( x )
        tmp = tmp[ !(tmp %in% levels) ]
        levels = c( levels, tmp )

      }

    }

  }

  type_check = F # Check if a correct coding scheme was provided
  if ( k > 1 ) {

    # Dummy coding
    if ( type == 'Dummy' | type == 'dummy' |
         type == 'D' | type == 'd' ) {

      # Initialize matrix
      m = matrix( 0, n, k - 1 )

      # Fill in values
      for ( l in 2:k ) {
        m[ x == levels[l], l - 1 ] = 1;
      }

      # Label columns
      colnames( m ) = paste( nm, as.character( levels[-1] ), sep = '_' )

      type_check = T
    }

    # Simple coding
    if ( type == 'Simple' | type == 'simple' |
         type == 'S' | type == 's' ) {

      # Initialize matrix
      m = matrix( -1/k, n, k - 1 )

      # Fill in values
      for ( l in 2:k ) {
        m[ x == levels[l], l - 1 ] = (k-1)/k;
      }

      # Label columns
      colnames( m ) = paste( nm, as.character( levels[-1] ), sep = '_' )

      type_check = T
    }

    # Effects coding
    if ( type == 'Effects' | type == 'effects' |
         type == 'Effect' | type == 'effect' |
         type == 'E' | type == 'e' ) {

      # Initialize matrix
      m = matrix( 0, n, k - 1 )

      # Fill in values
      for ( l in 2:k ) {
        m[ x == levels[l], l - 1 ] = 1;
      }
      m[ x == levels[1], ] = -1

      # Label columns
      colnames( m ) =  paste( nm, as.character( levels[-1] ), sep = '_' )

      type_check = T
    }

    # Intercept
    if ( type == 'Intercept' | type == 'intercept' |
         type == 'I' | type == 'i' ) {

      # Initialize matrix
      m = matrix( 0, n, k )

      # Fill in values
      for ( l in 1:k ) {
        m[ x == levels[l], l ] = 1;
      }

      # Label columns
      colnames( m ) =  paste( nm, as.character( levels ), sep = '_' )

      type_check = T
    }

    # Coef
    if ( type == 'Coefficient' |
         type == 'coefficient' |
         type == 'Coef' | type == 'coef' |
         type == 'C' | type == 'c' ) {

      # Check for custom values
      if ( is.null( weights ) ) {
        stop( 'Weights must be provided',
              call. = FALSE )
      }

      # Initialize matrix
      m = matrix( 0, n, 1 )

      # Fill in values
      for ( l in 1:k ) {
        m[ x == levels[l], 1 ] = weights[l];
      }

      # Label column
      colnames( m ) =  nm

      type_check = T
    }

  } else {
    # Intercept term
    m = matrix( 1, n, 1 );
    warning( 'Variable has only one level',
             call. = FALSE )
    type_check = T
  }

  if ( !type_check )
    stop( "Need type to be either 'Dummy', 'Simple', 'Effects',
          'Intercept', or 'Coef'.",
          call. = FALSE )

  return( m );
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
#' A function that calculates either Akaike's Information Criterion
#' (AIC; with or without a correction for small sample sizes) or the
#' Bayesian Information Criterion (BIC).
#'
#' @param logLik a log-likelihood value.
#' @param k the number of free parameters for the model.
#' @param n the number of observations in the sample.
#' @param type indicates whether to compute the 'AIC', 'AICc', or 'BIC'.
#'
#' @details
#'
#' Given a summed log-likelihood \eqn{L} from a model and \eqn{K} free
#' parameters, the AIC is \deqn{ 2K - 2L. }
#'
#' The AIC estimates the information loss from approximating the
#' true (generating) probability distribution with another probability
#' distribution. This discrepancy is represented by the Kullback-Leibler
#' information quantity, the negative of the generalized entropy.
#' Picking a model with the lowest information loss is asymptotically
#' equivalent to picking the model with the lowest AIC. Therefore,
#' the AIC is valid only for sufficiently large data sets.
#'
#' A correction for finite samples is recommended (e.g., Burnham &
#' Anderson, 2002), and for \eqn{N} observations the new equation is
#' \deqn{ 2K - 2L + \frac{2K(K+1)}{N+K+1}. } The corrected AIC is
#' denoted as 'AICc'.
#'
#' When comparing a set of candidate models, the AIC indicates which
#' model is most likely to best fit a new set of data generated from
#' the same process that produced the original sample. It is not
#' necessary to assume that the true generating model is part of the
#' set of candidate models. The AIC is subject to sampling variability;
#' a new sample of data will result in a different AIC value.
#'
#' The AIC has been criticized for being too liberal and likely to
#' select overly complex models. The AIC also neglects the sampling
#' variability of the parameter values. This means that if the
#' likelihood values for the parameters are not concentrated around
#' the maximum likelihood, using the AIC can lead to overly
#' optimistic assessments. The AIC is not consistent; as the number
#' of observations grows large, the probability of selecting a
#' true low-dimensional model based on model selection using the
#' AIC does not approach 1.
#'
#' The formula for the BIC is \deqn{ log(N)K - 2L. }
#'
#' The BIC is an asymptotic approximation to a Bayesian model
#' selection analysis. In Bayesian model selection, one must
#' compute the probability of each model given the data, which
#' requires specifying prior probabilities and integrating over
#' the parameter space. The BIC, though only an approximation,
#' is much easier to calculate and requires no specification of
#' priors. The BIC is consistent as the number of observations
#' grow large; the probability of selecting the true low-dimensional
#' model when using the BIC for model selection approaches 1.
#' The BIC also takes in account parameter uncertainty. However,
#' when using the BIC for model selection, one must assume that the
#' true generating model is in the set of candidate models
#' (which does not necessarily hold in reality).
#'
#' @references
#' Akaike, H. (1973). Information theory and an extension of the maximum
#'   likelihood principle. In B. N. Petrov & F. Caski (Eds.), Proceedings
#'   of the Second International Symposium on Information Theory (pp.
#'   267-281). Budapest:Akademiai Kiado.
#'
#' Burnham, K. P., & Anderson, D. R. (2002). Model selection and multimodel
#'   inference: A practical information-theoretic approach.
#'   New York:Springer-Verlag.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. Annals of
#'   Statistics, 6, 461-464.
#'
#' @return A value for either the AIC, AICc, or the BIC.
#'
#' @examples
#' N = 100; K = 2
#' x = rnorm( 100, 1, 2 )
#' m1 = sum( dnorm( x, 0, 1, log = T ) )
#' m2 = sum( dnorm( x, 1, 2, log = T ) )
#' # Corrected AIC values comparing the two models
#' print( round( infoCrit( c(m1,m2), K, N ), 2 ) )
#' # AIC values comparing the two models
#' print( round( infoCrit( c(m1,m2), K, N, type = 'AIC' ), 2 ) )
#' # BIC values comparing the two models
#' print( round( infoCrit( c(m1,m2), K, N, type = 'BIC' ), 2 ) )
#'
#' @export

infoCrit = function( logLik, k, n, type = 'AICc' ) {

  # Initialize output
  out = NA

  if ( type == 'AIC' ) out = AICc_val( logLik, k, n )
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

  string = paste(
    "  # Purpose:", "\n",
    "  # ...", "\n",
    "  # Arguments:", "\n",
    "  # ...", "\n",
    "  # Returns:", "\n",
    "  # ...", "\n", sep = '' )
  message( string )

}

# Lookup - 22
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

# Lookup - 26
#' Relative Weights for Information Criterion Values
#'
#' Computes the relative probability of being the
#' 'best' model for a set of candidate models based
#' on their AIC or BIC values.
#'
#' @param ic a vector of information criterion values,
#'   either a set of AIC or BIC values.
#'
#' @details
#'
#' A set of information criterion values (i.e., AIC or BIC)
#' can be transformed into a set of relative probabilities,
#' the probability that a given model is the 'best'. In the
#' case of the AIC, probabilities represent how likeliy
#' a given model will minimize information loss based on the
#' data and set of candidate models. In the case of the BIC,
#' probabilities represent how likely a given model is the
#' true (generating) model out of the set of candidate models
#' (assuming the true model is in the set).
#'
#' @return A vector of probabilities.
#'
#' @examples
#' # Model comparison with polynomial regression
#' x = rnorm(100) # Simulate predictor
#' df = data.frame( x1 = x, x2 = x^2, x3 = x^3, x4 = x^4 )
#' # True model is quadratic
#' df$y = .5 * df$x1 - .7 * df$x2 + rnorm(100,0,.75)
#' # Fit 4 models
#' m1 = lm( y ~ x1, data = df )
#' m2 = lm( y ~ x1 + x2, data = df ) # True
#' m3 = lm( y ~ x1 + x2 + x3, data = df )
#' m4 = lm( y ~ x1 + x2 + x3 + x4, data = df )
#' # Compute AIC
#' aic = c( AIC( m1 ), AIC( m2 ), AIC( m3 ), AIC( m4 ) )
#' names( aic ) = c( 'M1', 'M2', 'M3', 'M4' )
#' # AIC weights
#' print( round( icWeights( aic ), 2 ) )
#' # Compute BIC
#' bic = c( BIC( m1 ), BIC( m2 ), BIC( m3 ), BIC( m4 ) )
#' names( bic ) = c( 'M1', 'M2', 'M3', 'M4' )
#' # BIC weights
#' print( round( icWeights( bic ), 2 ) )
#'
#' # xa = seq( min(x), max(x), length = 100 )
#' # nd = data.frame( x1 = xa, x2 = xa^2, x3 = xa^3, x4 = xa^4, y = NA )
#' # plot( df$x1, df$y, pch = 19, xlab = 'x', ylab = 'y' )
#' # lines( xa, predict( m1, newdata = nd ), col = 'blue' )
#' # lines( xa, predict( m2, newdata = nd ), col = 'red' )
#' # lines( xa, predict( m3, newdata = nd ), col = 'green' )
#' # lines( xa, predict( m4, newdata = nd ), col = 'orange' )
#'
#' @export
icWeights = function( ic ) {

  delta_ic = ic - min( ic )
  relative_likelihood = exp( -.5 * delta_ic )
  weights = relative_likelihood / sum( relative_likelihood )

  return( weights )
}

# Lookup - 27
#' Simple Bootstrap Function
#'
#' A function for computing a test statistic over
#' multiple replicates of a set of data. Replicates
#' are created by sampling from the original data with
#' replacement.
#'
#' @param x a numeric vector, or a long-form matrix or data frame.
#' @param T_x a function to compute a desired test statistic
#'   (default is the mean).
#' @param nRep the number or replicates to sample.
#' @param ... additional parameters for the test statistic function.
#'
#' @return Returns a list containing the test statistics for each
#'   replicate, along with the test statistic for the original
#'   set of data.
#'
#' @examples
#' # Vector
#' x = rnorm( 100 )
#' bootstrap = simpleBootstrap( x, T_x = median )
#' quantile( bootstrap$replicates, c( .16, .84 ) )
#'
#' # Matrix
#' x = cbind( rlnorm( 100, -2, .5 ), rbeta( 100, 4, 1 ) )
#' bootstrap = simpleBootstrap( x, T_x = colMeans )
#' apply( bootstrap$replicates, 1, quantile, prob = c( .16, .84 ) )
#'
#' @export

simpleBootstrap = function( x, T_x = mean, nRep = 1000, ... ) {

  if ( is.vector(x) ) {
    # If x is a vector

    f = function( nr, dat, ... ) {

      # Permute data
      ord = 1:length(dat)
      y = x[ sample( ord, size = max(ord), replace = T ) ]

      # Compute test statistic
      out = T_x( y, ... )

      return( out )
    }

  } else if ( is.matrix(x) | is.data.frame(x) ) {
    # If x is a matrix or data frame

    f = function( nr, dat, ... ) {

      # Permute data
      ord = 1:dim(dat)[1]
      y = x[ sample( ord, size = max(ord), replace = T ), ]

      # Compute test statistic
      out = T_x( y, ... )

      return( out )
    }

  } else {
    # Return an error message

    string = 'Input should be a vector, or a long-form matrix or data frame.'
    stop( string, call. = FALSE )
  }

  #
  out = sapply( 1:nRep, f, dat = x, ... )

  # Return output
  return( list( replicates = out, observed = T_x( x, ... ) ) )
}

# Lookup - 28
#' Template for Project Headers
#'
#' Prints a template for creating a project header with the
#' author's name, email contact, and the date.
#'
#' @return Prints the template to the console window, where it
#'   can be copied and modified by the user.
#'
#' @examples
#' headerTemplate()
#'
#' @export

headerTemplate = function() {

  string = paste(
    '# Title', '\n',
    '# Written by Kevin Potter', '\n',
    '# email: kevin.w.potter@gmail.com', '\n',
    '# Please email me directly if you ', '\n',
    '# have any questions or comments', '\n',
    '# Last updated ', Sys.Date(),
    sep = '' )
  message( string )

}

# Lookup - 29
#' Add Horizontal Lines to a Plot
#'
#' Draws a set of horizontal lines all of the same
#' width onto an existing plot.
#'
#' @param yval a vector of y-axis values at which to
#'   draw the lines.
#' @param xl the left and right boundaries for the width of
#'   the lines.
#' @param ... additional parameters for the
#'   \code{\link[graphics]{segments}} function.
#'
#' @examples
#' plot( rnorm(100), rnorm(100), xlim = c(-3,3) )
#' horizLines( -2:2, c(-4,4) )
#'
#' @export

horizLines = function( yval, xl, ... ) {

  l = length( yval )

  segments( rep( xl[1], l ), yval,
            rep( xl[2], l ), yval, ... )

}

# Lookup - 30
#' Add Vertical Lines to a Plot
#'
#' Draws a set of vertical lines all of the same
#' height onto an existing plot.
#'
#' @param xval a vector of x-axis values at which to
#'   draw the lines.
#' @param yl the bottom and top boundaries for the height of
#'   the lines.
#' @param ... additional parameters for the
#'   \code{\link[graphics]{segments}} function.
#'
#' @examples
#' plot( rnorm(100), rnorm(100), ylim = c(-3,3) )
#' vertLines( -2:2, c(-4,4) )
#'
#' @export

vertLines = function( xval, yl, ... ) {

  l = length( xval )

  segments( xval, rep( yl[1], l ),
            xval, rep( yl[2], l ), ... )

}

# Lookup - 31
#' Fit a Basic Piecewise Regresssion
#'
#' Given an indepdent and dependent variable, along
#' with a vector of break points, fits a basic
#' piecewise regression to the data.
#'
#' @param x the independent (predictor) variable.
#' @param y the dependent variable.
#' @param breaks a vector of break points (in the same
#'   scale as the independent variable).
#' @param ... additional parameters for the
#'   \code{\link[stats]{lm}} function.
#'
#' @return A list with the output from the
#'   \code{\link[stats]{lm}} function,
#'   a data frame with the dummy coding used to
#'   construct the  piecewise regression, and
#'   x and y-axis values that can be used for
#'   plotting the estimated line segments.
#'
#' @examples
#' # Define generating parameters
#' beta = c( 1, 1, 1, -1 )
#' sigma = .5
#'
#' # Define break point
#' breaks = 0
#'
#' # Simulate data
#' x = runif( 100, -1, 1 ) # Predictor
#' # Define design matrix
#' X = matrix( 1, 100, 2 )
#' X[ x > breaks, 1] = 0; X[ x <= breaks, 2] = 0
#' X = cbind( X, X*x )
#' colnames( X ) = c( 'I1', 'I2', 'S1', 'S2' )
#' # Generate dependent variable
#' y = as.vector( X %*% cbind( beta ) )
#' y = y + rnorm( 100, 0, sigma )
#'
#' # Fit piecewise regression
#' fit = simplePiecewiseReg( x, y, breaks )
#' summary( fit$lm )
#'
#' Plot results
#' plot( x, y, pch = 19 )
#' lines( fit$xa, fit$ya )
#'
#' @export

simplePiecewiseReg = function( x, y, breaks, ... ) {

  # Number of breaks
  nb = length( breaks )
  nbp = nb + 1

  # Number of observations
  n = length( y )

  # Initialize matrix for x and y values
  M = matrix( NA, n, 2 + nbp*2 )
  M[,1] = x; M[,2] = y;

  elbows = data.frame(
    lb = c( min(x) - 1, breaks ),
    ub = c( breaks, max(x) ) )

  # Separate intercepts for each segment
  M[ , 1:nbp + 2 ] = apply( elbows, 1,
                            function(x)
                              as.numeric(
                                M[,1] > x[1] & M[,1] <= x[2] ) )
  # Separate slopes for each segment
  M[ , 1:nbp + nbp + 2 ] =
    M[ , 1:nbp + 2 ] * x

  # Define variable labels
  cn = c(
    'x', 'y',
    paste( 'I', 1:nbp, sep = '' ),
    paste( 'S', 1:nbp, sep = '' ) )
  colnames( M ) = cn

  # Convert to data frame
  M = as.data.frame( M )

  # Create formula for regression
  piecewise_reg_formula =
    paste( 'y ~ -1 + ',
           paste( cn[ grep( 'I', cn ) ], collapse = ' + ' ),
           ' + ',
           paste( cn[ grep( 'S', cn ) ], collapse = ' + ' ),
           collapse = '' )
  piecewise_reg_formula =
    as.formula( piecewise_reg_formula )

  # Fit piecewise regression
  fit = lm( piecewise_reg_formula, data = M, ... )

  # Output
  out = list(
    lm = lm( piecewise_reg_formula, data = M, ... ),
    data = M )

  # To plot unbroken segments, subsequent
  # segments following the first must
  # start from the previous value

  # Create indices to insert previous values
  # into data matrix
  break_pos = sapply( breaks, function(x)
    min( which( M[,1] >= x ) ) )
  break_pos = c( 1, break_pos, n )
  index1 = numeric( n + nb )
  index2 = numeric( n + nb )
  inc = 0
  for ( i in 2:length( break_pos ) ) {
    val2 = break_pos[i-1]:break_pos[i]
    val = val2
    if ( i > 2 ) val[1] = val2[1] + 1
    index1[1:length(val) + inc] = val
    index2[1:length(val) + inc] = val2
    inc = inc + length(val)
  }
  # Create new matrix of data to use to
  # generate x and y values for plotting
  # purposes
  xa = seq( range(x)[1], range(x)[2], length = 100 )

  # Initialize matrix for x and y values
  PM = matrix( NA, 100, nbp*2 )

  elbows = data.frame(
    lb = c( min(xa) - 1, breaks ),
    ub = c( breaks, max(xa) ) )

  # Separate intercepts for each segment
  PM[ , 1:nbp ] = apply( elbows, 1,
                         function(x)
                           as.numeric(
                             xa > x[1] & xa <= x[2] ) )
  # Separate slopes for each segment
  PM[ , 1:nbp + nbp ] =
    PM[ , 1:nbp ] * xa

  # Generate plotting values
  out$xa = xa
  out$ya = PM %*% cbind( coef( fit ) )

  return( out )
}

# Lookup - 32
#' Add Basic Color Map to a Plot
#'
#' @param Z a matrix with the numeric values
#'   to map to a set of colors, or a color
#'   map list object of class 'CM' (the
#'   output for this function).
#' @param vr an optional vector giving the lower and
#'   upper boundaries for to use when mapping the
#'   values to colors (defaults to the range for Z).
#' @param nc the number of equally-spaced
#'   intervals to use when mapping the
#'   values to colors (default approach needs an odd
#'   number).
#' @param clr the color values to which values should
#'   be mapped.
#' @param z the endpoints for the intervals to use
#'   when mapping the values to colors (need to be
#'   sequential and equal to the number of colors plus one).
#' @param border controls the border for each color square.
#' @param add logical; if true, adds the color map to an
#'   existing plot.
#' @param return logical; if true, returns a color map list
#'   object.
#' @param ... Additional inputs for the
#'   \code{\link[graphics]{polygon}} function.
#'
#' @return If indicated, returns a color map list object
#'   consisting of...
#'   \itemize{
#'     \item Z: the original matrix of data.
#'     \item cmap: a list with z, the endpoints
#'       of the intervals used to map colors to values,
#'       and clr, the set of colors to which values
#'       were assigned.
#'     \item x: a matrix with the x-axis values for
#'       each color cell for the \code{polygon} function.
#'     \item y: a matrix with the y-axis values for
#'       each color cell for the \code{polygon} function.
#'     \item clr: a list giving 1) the color for each cell
#'       in the map, 2) the x-axis position, 3) the y-axis
#'       position, and 4) the associated value.
#'     \item xl: the lower and upper intervals for the
#'       width of the color map.
#'     \item yl; the lower and upper intervals for
#'       the height of the color map.
#'   }
#'
#' @examples
#' # Generate a correlation matrix
#' n = 5 # Number of variables
#' A = matrix( runif(n^2)*2-1, ncol = n )
#' # Covariance matrix
#' S = t(A) %*% A
#' Z = cov2cor( S )
#'
#' # Create a blank plot
#' plot( c(0,5), c(0,5), type = 'n', xlab = 'X', ylab = 'Y' )
#' CM = addColorMap( Z )
#'
#' # Custom color map
#' plot( c(0,5), c(0,5), type = 'n', xlab = 'X', ylab = 'Y' )
#' addColorMap( Z, clr = heat.colors(11), return = F )
#'
#' # Custom intervals
#' plot( c(0,5), c(0,5), type = 'n', xlab = 'X', ylab = 'Y' )
#' z = c( -1, -.75, -.5, 0, .5, .75, 1 ); nc = length(z)
#' addColorMap( Z, z = z, nc = nc, return = F )
#'
#' # Create a color map initially, then add to plot
#' Z = matrix( runif(40), 10, 4 )
#' CM = addColorMap( Z, clr = heat.colors(11), add = F )
#' plot( CM )
#'
#' @export

addColorMap = function( Z,
                        vr = NULL,
                        nc = 11,
                        clr = NULL,
                        z = NULL,
                        border = NA,
                        add = T,
                        return = T,
                        ... ) {

  # If an existing color map list object
  # is provided
  if ( is.CM( Z ) ) {
    CM = Z

    # Dimensions
    nx = ncol( CM$Z )
    ny = nrow( CM$Z )
  }

  # If an existing color map list is
  # not provided
  if ( !is.CM( Z ) ) {

    # Dimensions
    nx = ncol( Z )
    ny = nrow( Z )

    # Range
    if ( is.null( vr ) ) {
      vr[1] = min( Z )
      vr[2] = max( Z )
    }

    # Default colors
    if ( is.null( clr ) ) {

      # Set to odd number
      if ( nc %% 2 != 1 ) {
        nc = nc + 1
        string = paste(
          'Default colors require odd number' )
        warning( string, call. = FALSE )
      }

      clr = c(
        rgb( seq( 0, .8, length = (nc-1)/2 ), 1, 1 ),
        rgb( 1, 1, 1 ),
        rev( rgb( 1, seq( 0, .8, length = (nc-1)/2 ), 1 ) ) )

    }

    # Default intervals to use for mapping values
    # to colors
    if ( is.null( z ) ) {

      z = seq( vr[1], vr[2], length = nc + 1 )

    }

    # Check that number of colors/intervals
    # are correct
    if ( !is.null( clr ) & !is.null( z ) ) {

      if ( !( length(clr) = length( z ) - 1 ) ) {

        string = paste(
          'Incorrect number of intervals to map',
          'values to colors - using defaults.' )
        warning( string, call. = FALSE )

        # Set to odd number
        if ( nc %% 2 != 1 ) {
          nc = nc + 1
          string = paste(
            'Default colors require odd number' )
          warning( string, call. = FALSE )
        }

        # Use default valus
        clr = c(
          rgb( seq( 0, .8, length = (nc-1)/2 ), 1, 1 ),
          rgb( 1, 1, 1 ),
          rev( rgb( 1, seq( 0, .8, length = (nc-1)/2 ), 1 ) ) )
        z = seq( vr[1], vr[2], length = nc + 1 )

      } else {
        nc = length( clr )
      }
    }

    # List of plotting variables
    CM = list(
      Z = Z,
      cmap = list(
        z = z,
        clr = clr ),
      x = cbind(
        rep( 0:(nx-1), each = ny ),
        rep( 0:(nx-1), each = ny ),
        rep( 1:nx, each = ny ),
        rep( 1:nx, each = ny ) ),
      y = cbind(
        rep( 0:(ny-1), nx ),
        rep( 1:ny, nx ),
        rep( 1:ny, nx ),
        rep( 0:(ny-1), nx ) )
    )

    # Specify color gradient
    CM$clr = data.frame( clr = rep( rgb( 0,0,0 ), nx * ny ),
                         x = NA, y = NA,
                         z = NA,
                         stringsAsFactors = F )

    # Loop over matrix cells
    inc = 1
    for ( i in 1:nx ) {
      for ( j in ny:1 ) {
        sel = sum( CM$Z[j,i] >= CM$cmap$z[ -( nc + 1 ) ] )
        CM$clr$clr[inc] =
          CM$cmap$clr[sel]
        CM$clr$x[inc] = i; CM$clr$y[inc] = j
        CM$clr$z[inc] = CM$cmap$z[sel]
        inc = inc + 1
      }
    }

    xl = c( 0, nx ); yl = c( 0, ny )
    CM$xl = xl; CM$yl = yl

    # Create a 'CM' class
    class( CM ) = 'CM'

  }

  # Add color map
  if ( add ) {

    # Loop over each cell
    for ( i in 1:(nx*ny) ) {

      polygon( CM$x[i,],
               CM$y[i,],
               col = CM$clr$clr[i],
               border = border,
               ... )

    }
  }

  # If indicated, return color map list object
  if ( return & !is.CM( Z ) ) return( CM )
}


#' @rdname addColorMap
#' @export
is.CM = function(x) inherits(x, "CM")

#' @rdname addColorMap
#' @export
print.CM = function( x, digits = 2, ... ) {

  nc = length( x$cmap$clr )
  out = data.frame(
    Colors = x$cmap$clr,
    Midpoints = x$cmap$z[-(nc+1)] + diff(x$cmap$z)/2,
    Lower = x$cmap$z[-(nc+1)],
    Upper = x$cmap$z[-1]
  )
  print( out, digits = digits, ... )

}

#' @rdname addColorMap
#' @export
plot.CM = function( x,
                    type = 'n',
                    xlab = ' ',
                    ylab = ' ',
                    xaxt = 'n',
                    yaxt = 'n',
                    bty = 'n',
                    ... ) {

  xl = x$xl; yl = x$yl
  plot( xl, yl,
        type = type,
        xlab = xlab,
        ylab = ylab,
        xaxt = xaxt,
        yaxt = yaxt,
        bty = bty,
        ... )
  addColorMap( x, return = F )

}

# Lookup - 33
#' Raise a Value to a Power
#'
#' This function raises a value x to a power a.
#' It is useful for notational purposes,
#' clearly separting the value and power from
#' other aspects of an equation.
#'
#' @param x a numerical value.
#' @param a the exponent.
#'
#' @return Returns the result of raising x to
#'   the power a.
#'
#' @examples
#' x = 2^4
#' y = pow( 2, 4 )
#' x == y
#'
#' # Sometimes it can be hard to determine
#' # what is being raised to a power
#' x = 2 * 13 ^ 4 / 8
#' # This function makes it easier to tell
#' x = 2 * pow( 13, 4 ) / 8
#'
#' @export

pow = function( x, a ) {

  out = x^a

  return( out )
}

# Lookup - 34
#' Custom Function to Standardize a Variable
#'
#' This function centers and scales a variable,
#' and adds attributes that allow the function
#' to convert the standardized variable back
#' into its original scale.
#'
#' @param x a numerical value.
#' @param reverse logical; if true, the function
#'   extracts the 'mean' and 'scale' attributes
#'   and reverses the standardization.
#'
#' @return If reverse is FALSE, returns the standardized
#' values for x, with attributes giving the mean and scale;
#' otherwise, returns the raw, unstandardized values
#' of x.
#'
#' @examples
#' x = runif( 30 )
#' # Standardize the variable
#' z = my_standardize( x )
#' # Unstandardize the variable
#' y = my_standardize( x, reverse = T )
#' all( x == y )
#'
#' @export

my_standardize = function( x, reverse = F ) {

  if ( !reverse ) {
    m = mean( x, na.rm = T )
    s = sd( x, na.rm = T )
    out = ( x - m ) / s
    attributes( out ) = list( mean = m, scale = s )
  } else {
    m = attributes( x )$mean
    s = attributes( x )$scale
    out = x * s + m
    attributes( out ) = NULL
  }

  return( out )
}

# Lookup - 35
#' Functions to Evaluate Run Times
#'
#' The functions \code{tick} and \code{tock} can be
#' used to roughly estimate the run time of a segment
#' of code.
#'
#' @return The function \code{tick} creates a global
#' variable 'run_time' with the current system time.
#' The function \code{tock} then updates 'run_time'
#' with the difference between the current system
#' time and the previous value of 'run_time'.
#'
#' @examples
#' tick()
#' Sys.sleep(2.5)
#' tock()
#' print( run_time )
#'
#' @export

tick = function() {

  # Create a global variable with
  # the current time
  run_time <<- Sys.time()

}

#' @rdname tick
#' @export
tock = function() {

  if ( exists( 'run_time' ) ) {
    run_time <<- Sys.time() - run_time
  } else {
    warning( paste(
      "Variable 'run_time' does not exist;",
      "Use 'tick()' to create variable." ))
  }

}

# Lookup - 36
#' Function to Identify NA Values by Row
#'
#' This function identifies rows in a matrix
#' or data frame that contain any NA values.
#'
#' @param x a matrix or data frame.
#' @param any logical; if \code{TRUE} check
#'   if any observation in the row is an \code{NA}
#'   value, else checks if all observations are
#'   \code{NA}.
#'
#' @return A logical vector with values of
#' \code{TRUE} for rows with any \code{NA}
#' values (or rows where all values are \code{NA}
#' if \code{any} is \code{FALSE}).
#'
#' @examples
#' x = matrix( rnorm(9), 3, 3 )
#' x[2,3] = NA
#' findNA( x )
#' x = data.frame( A = c( 1, 2, NA ), B = 0 )
#' findNA( x )
#'
#' #' x = matrix( rnorm(9), 3, 3 )
#' x[2,] = NA
#' x[3,1] = NA
#' findNA( x, any = F )
#'
#' @export

findNA = function( x, any = T ) {

  # Initialize output
  out = NULL

  # If input is matrix or data frame
  if ( is.matrix( x ) |
       is.data.frame( x ) ) {
    # Check whether NA values are present in
    # any of the rows
    if ( any ) {
      out = apply( x, 1, function(r) any( is.na(r) ) )
    } else {
      out = apply( x, 1, function(r) all( is.na(r) ) )
    }
  } else {
    # Throw an error message
    string = "Input should be a matrix or data frame"
    stop( string, call. = F )
  }

  return( out )
}

# Lookup - 37
#' Transform Hit/False Alarm Rates into SDT Parameters
#'
#' Calculates d' and c parameter estimates for the
#' equal-variance Signal Detection Theory (SDT) model
#' for binary data using the algebraic method.
#' Also computes alternative metrics, including
#' A', B, and B'' (Stanislaw & Todorov, 1993).
#'
#' @param dat either 1) a vector of four values, the
#'   frequencies for hits and false alarms followed
#'   by the associated total number of trials for each,
#'   or 2) a vector of two values, the proportion of
#'   hits and false alarms.
#' @param centered logical; if \code{TRUE} uses the
#'   parameterization in which the distributions are
#'   equidistant from 0.
#' @param correct the type of correction to use, where...
#'   \itemize{
#'     \item 0 = none;
#'     \item 1 = The log-linear approach, adds .5 to
#'     the hits andfalse alarm frequencies, then adds
#'     1 to the total number of trials (Hautus, 1995);
#'     \item 2 = The conditional approach, where only
#'     proportions equal to 0 or 1 are adjusted by
#'     .5/N or (N-.5)/N respectively, where N is the
#'     associated number of total trials for the given
#'     proportion (Macmillan & Kaplan, 1985).
#'   }
#'
#' @return A named vector with five values: 1) d', the estimate
#' of separation between the noise and signal distributions;
#' 2) c, the estimate of response bias (the cut-off determining
#' whether a response is 'Signal' versus 'Noise'); 3) A', a
#' non-parametric estimate of discrimination; 4) B, the ratio
#' of whether a person favors responding 'Signal' over
#' whether he or she favors responding 'Noise'; 5) B'', a
#' non-parametric estimate of B.
#'
#' @references
#' Hautus, M. J. (1995). Corrections for extreme proportions and their
#'   biasing effects on estimated values of d'. Behavior Research Methods
#'   Instruments, & Computers, 27(1), 46 - 51. DOI: 10.3758/BF03203619
#'
#' Macmillan, N. A. & Kaplan, H. L. (1985). Detection theory analysis
#'   of group data: estimating sensitivity from average hit and
#'   false-alarm rates. Psychological Bulletin, 98(1), 185 - 199.
#'
#' Stanislaw, H. & Todorov, N. (1993). Calculation of signal detection
#'   theory measures. Behavior Research Methods, Instruments, & Computers,
#'   31, 137 - 149.
#'
#' @examples
#' # Hit/false alarm rate when d' = 1 and c = 0
#' x = c( H = pnorm( 0, -.5 ), FA = pnorm( 0, .5 ) )
#' round( binarySDT( x ), 2 )
#' # Hit/false alarm rate when d' = 1 and c = .25
#' x = c( H = pnorm( .25, -.5 ), FA = pnorm( .25, .5 ) )
#' round( binarySDT( x ), 2 )
#' # Using frequencies
#' y = c( round( x*20 ), 20, 20 )
#' round( binarySDT( y ), 2 )
#' # Correction for extreme values
#' y = c( 10, 6, 10, 10 )
#' round( binarySDT( y, correct = 1 ), 2 )
#' @export

binarySDT = function( dat, centered = T, correct = 0 ) {

  # Initialize values
  out = c( "d'" = NA, "c" = NA, "A'" = NA, "B" = NA, "B''" = NA )

  H = NULL; FA = NULL

  # Extract data

  # If counts and trial numbers are provided
  if ( length( dat ) == 4 ) {

    # Extract frequencies
    S = dat[1]; N = dat[2];
    # Extract total number of trials
    Ns = dat[3]; Nn = dat[4]

    # Determine hits and false-alarm rates
    H = S/Ns; FA = N/Nn;

    # Apply corrections if indicated
    if ( correct == 1 ) {
      H = (S+.5)/(Ns+1)
      FA = (N+.5)/(Nn+1)
    }
    if ( correct == 2 ) {
      if ( H == 0 ) H = .5/Ns
      if ( FA == 0 ) FA = .5/Nn
      if ( H == 1 ) H = (Ns-.5)/Ns
      if ( FA == 1 ) FA = (Nn-.5)/Nn
    }

  }

  # If only proportions are provided
  if ( length( dat ) == 2 ) {

    # Extract hits and false-alarm rates
    H = dat[1]; FA = dat[2]

    if ( correct > 0 ) {
      string = paste(
        'Correction cannot be applied',
        'using proportions; frequencies',
        'and total number of trials are',
        'required.' )
      warning( string, call. = F )
    }
  }

  if ( !is.null(H) & !is.null(FA) ) {

    if ( !centered ) {

      # Obtain estimate of d'
      dp_est = qnorm( H ) - qnorm( FA )

      # Obtain estimate of criterion
      k_est = qnorm( 1 - FA )

    } else {

      # Obtain estimate of criterion
      k_est = -.5*( qnorm( H ) + qnorm( FA ) )

      # Obtain estimate of d'
      dp_est = 2*( qnorm( H ) + k_est )

    }

    # Compute additional values
    num = pow( H - FA, 2 ) + (H - FA)
    denom = 4*max(H,FA) - 4*H*FA
    Ap_est = sign( H - FA )*(num/denom)

    log_beta_est = .5*( pow( qnorm( FA ), 2 ) - pow( qnorm( H ), 2 ) )

    num = H*(1-H) - FA*(1-FA)
    denom = H*(1-H) + FA*(1-FA)
    beta_pp = sign( H - FA)*( num/denom )

    out[1] = dp_est
    out[2] = k_est
    out[3] = Ap_est
    out[4] = exp( log_beta_est )
    out[5] = beta_pp

  }

  return( out )
}

# Lookup - 38
#' Estimate Parameters for a Two-Boundary Unbiased Wiener Process
#'
#' Estimates the drift rate, boundary separation, and
#' non-decision time for a two-boundary wiener process
#' with an equidistant starting point using a R function
#' adapted from Wagenmaker et al. (2007).
#'
#' @param dat a vector of inputs:
#'   \enumerate{
#'     \item Total number of trials;
#'     \item Number of correct responses;
#'     \item Mean response time for correct responses;
#'     \item Variance for corrct response times;
#'   }
#' @param s the scaling parameter (typically set to either 0.1 or 1)
#'
#' @return An estimate of the drift rate, boundary separation,
#' and non-decision time. Bias is fixed to 0.5.
#'
#' @references
#' Wagenmakers, E. J., Van Der Maas, H. L., & Grasman, R. P.
#'   (2007). An EZ-diffusion model for response time and accuracy.
#'   Psychonomic bulletin & review, 14, 3-22.
#'   https://doi.org/10.3758/BF03194023.
#'
#' @examples
#' x = c( N = 100, Correct = 80, MRT = .8, VRT = .4^2 )
#' EZdiffusion( x )
#' @export

EZdiffusion = function( dat, s = 1 ) {

  # Initialize output
  out = c(
    drift_rate = NA,
    bias = NA,
    boundary_sep = NA,
    non_decision_time = NA
  )

  # Extract data
  Trials = dat[1]
  Correct = dat[2]
  MRT = dat[3]
  VRT = dat[4]

  # Check for invalid inputs
  invalid_inputs = c(
    # Missing data
    any( is.na( dat ) ),
    # Non-positive mean response time
    MRT <= 0,
    # Non-positive variance
    VRT <= 0,
    # Number of trials is not an integer
    Trials %% 1 != 0,
    # Number of correct responses is not an integer
    Correct %% 1 != 0,
    # Number of correct responses exceeds number of trials
    Correct > Trials
  )

  invalid_inputs_warning = c(
    'NA values',
    'Non-positive mean response time',
    'Non-positive variance for response time',
    'Non-integer value for trials',
    'Non-interger value for frequency correct',
    'Frequency correct exceeds number of trials'
  )

  # Return NA for invalid inputs
  if ( any( invalid_inputs ) ) {
    string = paste(
      'Invalid inputs:\n',
      paste( invalid_inputs_warning[invalid_inputs], collapse = '\n' ),
      sep = '' )
    warning( string, call. = F )
    return( out )
  }

  # Proportion of correct responses
  Pc = Correct / Trials

  # Scaling parameter
  s2 = pow( s, 2 )

  # No information if P( Correct ) = 0, 0.5, or 1
  if ( Pc %in% c( 0, .5, 1 ) ) {

    warning( paste(
      'Insufficient information for estimation',
      ' - an edge correction will be applied' ),
      call. = FALSE )

    # For the edge correction, adjust results
    # by one half of an error:
    if ( Pc == 1 )
      Pc = 1 - 1 / ( 2*Trials)
    if ( Pc == 0 )
      Pc = 1 + 1 / ( 2*Trials)
    if ( Pc == .5 )
      Pc = .5 + 1 / ( 2*Trials)

  }

  # Compute the log odds for P( Correct )
  L = qlogis(Pc)
  x = L * ( L * pow( Pc, 2 ) - L * Pc + Pc - .5 )/VRT

  # Estimate drift and boundary separation
  # from P( Correct ) and the variance of correct RT
  drift_rate = sign( Pc - .5 ) * s * pow( x, 1/4 )
  boundary_sep = s2 * L / drift_rate

  # Estimate the non-decision time
  y = -drift_rate * boundary_sep/s2
  MDT = ( boundary_sep/( 2 * drift_rate) ) *
    ( 1 - exp(y) )/( 1 + exp(y) )
  non_decision_time = MRT - MDT

  # Return estimates of parameters
  # if ( Correct / Trials == .5 ) drift_rate = 0
  out[1] = drift_rate
  out[2] = 0.5
  out[3] = boundary_sep
  out[4] = non_decision_time

  return( out )
}

# Lookup - 39
#' Plot a Heatmap for a Correlation Matrix
#'
#' Generates a heatmap of the upper triangle of a
#' correlation matrix.
#'
#' @param df a data frame with all variables to include
#'   in the correlation matrix.
#' @param ttl an optional title for the figure.
#' @param labels the labels for the rows/columns.
#' @param new logical; if \code{TRUE} generates a
#'   new plotting window.
#' @param lyt an optional matrix specifying the
#'   layout of the main panel (1) versus
#'   the side panel (2) with the color gradient.
#' @param gradient the final end colors for
#'   the negative and positive correlations, respectively.
#' @param txtSz the size of the text in the figure.
#' @param mc_adjust the method to use when correcting for
#'   multiple comparisons (see ).
#' @param cut_off cut-off for significance.
#' @param H the height to use if a new plotting window is
#'   generated.
#' @param W the width to use if a new plotting window is
#'   generated.
#'
#' @return A heatmap for the upper-triangle portion of the
#'   correlation matrix.
#'
#' @examples
#' # Load data
#' data("mtcars")
#' my_data <- mtcars[, c(1,3,4,5,6,7)]
#' # Generate heatmap
#' correlationHeatmap( my_data, labels = colnames(my_data) )
#'
#' @export

correlationHeatmap <- function( df,
                                ttl = "Correlation matrix",
                                labels = NULL,
                                new = T, lyt = NULL,
                                gradient = c("orange", "blue"),
                                txtSz = 1.25, mc_adjust = 'BH',
                                cut_off = 0.05, H = 20/3, W = 25/3 ) {

  # Compute correlations
  omega = cor( df )

  # Compute p-values for correlations
  p_mat = matrix( NA, nrow( omega ), ncol( omega ) )
  p_val = rep( NA, sum( upper.tri( omega ) ) )
  k = 1
  for ( i in 1:nrow( omega ) ) {
    for( j in 1:ncol( omega ) ) {
      if ( i != j ) {
        tst = cor.test(
          df[[ rownames( omega )[i] ]],
          df[[ colnames( omega )[j] ]]
        )
        p_mat[i,j] = tst$p.value
        if ( j > i ) {
          p_val[k] = tst$p.value
          k = k + 1
        }
      }
    }
  }

  # Adjust for multiple comparisons
  p_val = p.adjust( p_val, method = mc_adjust )
  k = 1
  for ( i in 1:nrow( omega ) ) {
    for ( j in 1:ncol( omega ) ) {
      if ( j > k ) {
        p_mat[i,j] = p_val[k]
        p_mat[j,i] = p_mat[i,j]
        k = k + 1
      }
    }
  }
  diag( p_mat ) = 0

  # Create new plotting window
  if (new) {
    x11( height = H, width = W )
  }

  # Create panels for plot and color map
  nr = nrow(omega)
  pos = seq(nrow(omega), 0, -1)
  if (is.null(lyt)) {
    lyt = cbind(1, 1, 1, 1, 2)
  }
  layout(lyt)

  # Create blank plot
  xl = range(pos)
  yl = range(pos)
  mrg = c(11, 11, 0.5, 0.5)
  par(mar = mrg)
  blankPlot(xl, yl)

  # Function for drawing colored box
  draw_box = function(ri, ci, clr = NULL, brd = NULL, slash = FALSE ) {
    xa = c(ci + 1, ci + 1, ci, ci)
    ya = c(ri, ri + 1, ri + 1, ri)
    if (is.null(clr))
      clr = "white"
    if (is.null(brd))
      brd = clr
    if ( !slash ) {
      polygon(rev(pos)[xa], pos[ya], col = clr, border = brd)
    } else {
      segments( rev(pos)[xa][1],
                pos[ya][1],
                rev(pos)[xa][3],
                pos[ya][2] )
    }
  }

  # Color gradient
  r_range = seq(0, 1, 0.01)
  color_f = colorRampPalette(c("white", gradient[2]))
  color_pos = color_f(length(r_range))
  color_f = colorRampPalette(c("white", gradient[1]))
  color_neg = color_f(length(r_range))

  # Loop over upper triangle of correlation matrix
  for (ri in 1:nr) {
    for (ci in 1:nr) {
      cur_R = round(omega[ri, ci], 2)
      if (cur_R > 0) {
        cur_clr = color_pos[round(r_range, 2) == round(cur_R, 2)]
      }
      if (cur_R < 0) {
        cur_clr = color_neg[round(r_range, 2) == round(abs(cur_R), 2)]
      }
      if ( cur_R == 0 ) {
        cur_clr = 'white'
      }
      if (ri == ci) {
        cur_clr = "white"
      }
      if (ri >= ci) {
        cur_clr = "white"
      }
      cur_brd = NULL
      draw_box(ri, ci, clr = cur_clr, brd = cur_brd)
    }
  }

  # Denote non-significant boxes
  if (!is.null(p_mat)) {

    for (ri in 1:nr) {
      for (ci in 1:nr) {

        cur_R = round(omega[ri, ci], 2)
        cur_p = p_mat[ri, ci]
        if (is.na(cur_p))
          cur_p = 0

        if (cur_p > cut_off) {
          if (ri < ci) {
            if (cur_R > 0) {
              cur_clr = color_pos[round(r_range, 2) ==
                                    round(cur_R, 2)]
            }
            if (cur_R < 0) {
              cur_clr = color_neg[round(r_range, 2) ==
                                    round(abs(cur_R), 2)]
            }

            cur_brd = "black"
            draw_box(ri, ci, clr = cur_clr, brd = cur_brd, slash = T )
          }
        }
      }
    }
    draw_box(nr, 1, clr = "white", brd = "black")
    draw_box(nr, 1, clr = "white", brd = "black", slash = T )
    text(1.25, 0.5, paste("Non-significant at p >", cut_off),
         pos = 4, cex = txtSz)
  }

  if ( is.null( labels ) ) {
    labels = colnames( omega )
  }


  for (ri in 1:nr) {
    text(nr - (ri - 1), ri - 0.5, rev(labels)[ri], pos = 2,
         xpd = NA, cex = txtSz)
  }
  mtext(ttl, side = 1, line = 1, cex = txtSz)
  lbl = lowerUpper(0.1, omega[lower.tri(omega)])
  lbl = seq(lbl[1], lbl[2], 0.1)
  lbl = rev(lbl)
  pos = rev(0:(length(lbl)))
  xl = c(0, 2)
  yl = range(pos)
  par(mar = c(mrg[1], 0, mrg[3], mrg[4]))
  blankPlot(xl, yl)
  for (ri in 1:length(lbl)) {
    cur_R = lbl[ri]
    if (cur_R > 0) {
      cur_clr = color_pos[round(r_range, 2) == round(cur_R,
                                                     2)]
    }
    if (cur_R < 0) {
      cur_clr = color_neg[round(r_range, 2) == round(abs(cur_R),
                                                     2)]
    }
    if (cur_R == 0) {
      cur_clr = "white"
    }
    draw_box(ri, 1, clr = cur_clr)
  }
  string = as.character(lbl)
  string[lbl > 0] = paste(" ", string[lbl > 0])
  string[lbl == 0] = " 0.0"
  text(rep(txtSz, length(lbl)), pos[-1] + 0.5, string)
  mtext("R", side = 1, line = 1, cex = txtSz)

}

# Lookup - 40
#' Template for Standard Plotting Code
#'
#' This function prints to the console a
#' a template for an R function to generate a
#' standard plot.
#'
#' @return Template of an R function to generate a standard plot.
#'
#' @export

plotTemplate = function() {

  string = "
  plot_f = function( df,
                     inc,
                     lbl = c( 'x-axis', 'y-axis' ),
                     lnSz = 2,
                     ptSz = 1.25,
                     axSz = 1.25,
                     axPos = -1.5,
                     lblSz = 1.15,
                     lblPos = 1.5,
                     new = T ) {
    # Purpose:
    # ...
    # Arguments:
    # df     - A data frame with the x and y values
    # inc    - The spacing value for the x and y-axis limits,
    #          respectively
    # lbl    - The labels for the x and y-axis, respectively
    # lnSz   - The width of the lines
    # ptSz   - The size of the points
    # axSz   - The text size for the axis values
    # axPos  - The position of the axis values
    # lblSz  - The text size for the axis labels
    # lblPos - The position of the axis labels
    # new    - Logical; if TRUE, a new plotting window is created

    # Setup
    x = df$x
    y = df$y

    # If specifed, create new plotting window
    if ( new ) x11()

    # x and y-axis boundaries
    xl = lowerUpper( inc[1], x )
    yl = lowerUpper( inc[2], y )

    # Position of axis labels
    xax = seq( xl[1], xl[2], inc[1] )
    yax = seq( yl[1], yl[2], inc[2] )

    # Create blank plot
    blankPlot( xl, yl )

    # Grid lines
    horizLines( yax, xl, col = 'grey', lwd = lnSz )

    # Plotting
    lines( x, y, lwd = lnSz )
    points( x, y, pch = 19, cex = ptSz )

    # Axes and labels
    customAxes( xl, yl, lnSz = lnSz )

    # x-axis
    axis( 1, xax,
          tick = F, line = axPos, cex.axis = axSz )
    mtext( lbl[1],
           side = 1, line = lblPos, cex = lblSz )

    # y-axis
    axis( 2, yax,
          tick = F, line = axPos, cex.axis = axSz )
    mtext( lbl[2],
           side = 2, line = lblPos, cex = lblSz )

  }"

  message( string )
}

# Lookup - 41
#' Print a Nicely Formatted Table
#'
#' Given a data frame, prints to console
#' a formated table of the results.
#'
#' @param tbl a formatted data frame.
#' @param return logical; if \code{TRUE}, returns
#'   the character vector with each formatted
#'   row.
#'
#' @return If \code{return} is \code{TRUE} returns the
#' vector of character strings with the formatted output
#' for each row of the table.
#'
#' @examples
#' data( 'mtcars' )
#' tbl = aggregate( mtcars$mpg, mtcars[,c('cyl','vs')], mean )
#' tbl$x = round( tbl$x, 1 )
#' colnames( tbl ) = c( "# of cylinders", "Engine type", "Miles per gallon" )
#' printTable( tbl )
#'
#' @export

printTable = function( tbl, return = F ) {

  # Initialize output
  out = matrix( " ", nrow( tbl ) + 1, ncol( tbl ) )

  # Loop over columns of table
  for ( i in 1:ncol( tbl ) ) {

    # Determine maximum number of characters for elements
    # in the table's column (including the column name)
    nc = max( c( sapply( as.character( tbl[[i]] ), nchar ),
                 nchar( colnames(tbl)[i] ) ) )

    # Loop over the table's rows
    for ( j in 1:( nrow( tbl ) + 1 ) ) {

      if ( j > 1 ) {
        # Elements in column
        val = as.character( tbl[[i]] )[j-1]
      } else {
        # Column name
        val = colnames( tbl )[i]
      }
      # Current number of characters
      cur_nc = nchar( val )
      # If necessary pad out characters with empty spaces
      if ( cur_nc < nc ) {
        val = paste( paste( rep( " ", nc - cur_nc ), collapse = "" ),
                     val, sep = "" )
        out[j,i] = val
      } else {
        out[j,i] = val
      }

    }
  }

  # Convert to vector of strings
  output = apply( out, 1, paste, collapse = " | " )
  output = sapply( output, function(x) paste( x, "|" ) )

  if ( !return ) {
    for ( i in 1:length( output ) ) {
      cat( c( output[i], '\n' ) )
    }
  } else {
    return( output )
  }

}

# Lookup - 42
#' Estimate or Format p-values from Monte Carlo Samples
#'
#' Given a set of Monte Carlo samples, estimates a p-value
#' from the proportion of values that fall above or below
#' a comparison point. If \code{string} is \code{TRUE},
#' takes a numeric p-value and converts it into a
#' formatted character string, either 'p = ...' or
#' 'p < ...'.
#'
#' @param x either 1) a vector of numeric values (Monte
#'   Carlo samples) or 2) a single p-value.
#' @param comparison the comparison point used to compute
#'   the p-value for the Monte Carlo samples.
#' @param alternative a character string, either 1) 'two-sided',
#'   2) 'greater', or 3) 'less', indicating the type of
#'   alternativ hypothesis to test. If 'two-sided', uses the
#'   alternativ hypothsis that produces the smallest p-value,
#'   and then adjusts for the two-sided comparison by multiplying
#'   by two.
#' @param digits the number of digits to round to when formatting
#'   the numeric p-value
#'
#' @return Either a numeric p-value or a character string,
#' a nicely formatted version of the p-value.
#'
#' @examples
#' x = rnorm( 1000 )
#' p = pvalueMC( x )
#' print( pvalueMC( p ) )
#' p = pvalueMC( x, alternative = 'greater' )
#' print( pvalueMC( p ) )
#' p = pvalueMC( x, alternative = 'less' )
#' print( pvalueMC( p, digits = 4 ) )
#'
#' @export

pvalueMC = function( x, comparison = 0,
                     alternative = 'two-sided',
                     digits = 3 ) {

  # Initialize output
  out = NULL

  # Assume if siz of x equals 1 that
  # it is a numeric p-value
  string = FALSE
  if ( length( x ) == 1 ) {

    if ( x >= 0 &
         x <= 1 ) {

      string = TRUE
    }
  }

  # Estimate p-value from Monte Carlo samples
  if ( !string ) {

    check = FALSE

    # Two-sided test
    if ( alternative == 'two-sided' ) {
      check = TRUE

      if ( median( x ) > comparison ) {
        out = mean( x < comparison )
      } else {
        out = mean( x > comparison )
      }
      out = out*2

    }

    # Test if greater than comparison
    if ( alternative == 'greater' ) {

      check = TRUE

      out = mean( x < comparison )

    }

    # Test if less than comparison
    if ( alternative == 'less' ) {

      check = TRUE

      out = mean( x > comparison )

    }

    # Informative error message if
    # alternativ misspecified
    if ( !check ) {

      err_msg = paste(
        "Please specify 'alternative' as",
        "either 'two-sided', 'greater', or",
        "'less'." )
      stop( err_msg )

    }

  } else {

    # Convert numeric p-value into
    # a nice string character
    p = round( x, digits )
    out = paste( "p = ", p, sep = "" )
    if ( p == 0 ) {
      nd = digits - 1
      nd = paste(
        "0.",
        paste( rep( 0, nd ), collapse = '' ),
        '1', sep = '' )
      out = paste( "p < ", nd, sep = "" )
    }

  }

  return( out )
}

# Lookup - 43
#' Colorblind Friendly Palette
#'
#' Generates a set of 8 colors that are
#' colorblind friendly.
#'
#' @references
#' See http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette.
#'
#' @return A vector of 8 hexadecimal strings giving
#' # colors (starting with grey) that are colorblind
#' friendly.
#'
#' @examples
#' clrs = colorblindPalette()
#' plot( 1:8, 1:8, pch = 15, cex = 4, col = clrs )
#'
#' @export

colorblindPalette = function() {

  # The palette with grey:
  out = c(
    grey = "#999999",
    orange = "#E69F00",
    light_blue = "#56B4E9",
    green = "#009E73",
    yellow = "#F0E442",
    blue = "#0072B2",
    red = "#D55E00",
    pink = "#CC79A7" )

  return( out )
}

