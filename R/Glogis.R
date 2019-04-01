#' @title Generalized Logistic Distribution
#' @description Family objects for generalized logistic distribution provide a convenient way to specify the details of the models used by functions such as \code{elliptical}.
#' @param parma parameter vector (alpha, m).
#' @return An object of class "family" for generalized logistic distribution.
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
#' elliptical.fitGl <- elliptical(y ~ x1+x2+x3, family = Glogis(c(1,1))
#' ,data=luz)
#' family(elliptical.fitGl)
#' @export


Glogis <- function (parma = stop("no alpha=alpha(m) or m argument")) 
{
  if ((parma[1] <= 0) || (parma[2] <= 0)) 
    stop(paste("alpha=alpha(m) and m must be positive"))
  make.family.elliptical("Glogis", arg = list(alpha = parma[1], 
                                              mp = parma[2]))
}
