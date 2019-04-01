#' @title Logistic Type I Distribution
#' @description Family objects for logistic type I distribution provide a convenient way to specify the details of the models used by functions such as \code{elliptical}.
#' @return An object of class "family" for logistic type I distribution.
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
#' elliptical.fitLI <- elliptical(y ~ x1+x2+x3, family = LogisI()
#' ,data=luz)
#' family(elliptical.fitLI)
#' @export


LogisI <- function () 
{
  make.family.elliptical("LogisI")
}