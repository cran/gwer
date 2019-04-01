#' @title Contaminated Normal Distribution
#' @description Family objects for contaminated normal distribution provide a convenient way to specify the details of the models used by functions such as \code{elliptical}.
#' @param parmt parameters vector (epsi, sigma).
#' @return An object of class "family" for contaminated normal distribution.
#' @references Fang, K. T., Kotz, S. and NG, K. W. (1990, ISBN:9781315897943).
#' Symmetric Multivariate and Related Distributions. London: Chapman and Hall.
#' @seealso \code{\link{family.elliptical}}, \code{\link{elliptical}}, \code{\link{gwer}}
#' @keywords elliptical distributions
#' @examples
#' \dontrun{
#' data(luzdat)
#' y <- luzdat$y
#' x1 <- luzdat$x1 ; x1 <- factor(x1) ; x1 <- C(x1,treatment)
#' x2 <- luzdat$x2
#' x3 <- (luzdat$x2)^2
#' luz <- data.frame(y,x1,x2,x3)
#' elliptical.fitCn <- elliptical(y ~ x1+x2+x3, family = Cnormal(c(1,.5))
#' ,data=luz)
#' family(elliptical.fitCn)
#' }
#' @export


Cnormal <- function (parmt = stop("no epsi or sigma argument")) 
{
  stop(paste("not implement yet"))
  if ((parmt[1] < 0) || (parmt[1] > 1) || (parmt[2] <= 0)) 
    stop(paste("0<=epsilon<=1 and sigma must be positive"))
  make.family.elliptical("Cnormal", arg = list(epsi = parmt[1], 
                                               sigmap = parmt[2]))
}
