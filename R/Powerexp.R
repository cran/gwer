#' @title Power Exponential Distribution
#' @description Family objects for power exponential distribution provide a convenient way to specify the details of the models used by functions such as \code{elliptical}.
#' @param k shape parameter.
#' @return An object of class "family" for power exponential distribution.
#' @references Fang, K. T., Kotz, S. and NG, K. W. (1990, ISBN:9781315897943).
#' Symmetric Multivariate and Related Distributions. London: Chapman and Hall.
#' @seealso \code{\link{family.elliptical}}, \code{\link{elliptical}}, \code{\link{gwer}}
#' @keywords elliptical distributions
#' @examples
#' data(luzdat)
#' y <- luzdat$y
#' x1 <- luzdat$x1 ; x1 <- factor(x1) ; x1 <- C(x1,treatment)
#' x2 <- luzdat$x2
#' x3 <- (luzdat$x2)^2
#' luz <- data.frame(y,x1,x2,x3)
#' elliptical.fitPe <- elliptical(y ~ x1+x2+x3, family = Powerexp(1)
#' ,data=luz)
#' family(elliptical.fitPe)
#' @export


Powerexp <- function (k = stop("no k argument")) 
{
  if (abs(k) > 1) 
    stop(paste("k must be (-1,1)"))
  make.family.elliptical("Powerexp", arg = k)
}
