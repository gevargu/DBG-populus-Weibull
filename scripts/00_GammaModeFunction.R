##' Jenkins
##'
##' This code calculates Shape parameters from central tendency and scale for a Gamma distribution using median and standard deviation.
##' 
##' @references Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: Tutorial with R, JAGS, and Stan. Academic Press, Elsevier.
##' 
##' @param mode vector with 
##' @param sd column name as character with the date identifier.
##' 
##' @author John K. Kruschke
##' @export
##' @examples
##' \dontrun{
##' #nothing yet
##' }
gammaShRaFromModeSD = function( mode , sd ) {
  if ( mode <=0 ) stop("mode must be > 0")
  if ( sd <=0 ) stop("sd must be > 0")
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}