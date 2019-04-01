#' @title Generalized Student Distribution
#' @description Family objects for generalized student distribution provide a convenient way to specify the details of the models used by functions such as \code{elliptical}.
#' @param parm parameter vector (s, r) for this distribuition.
#' @return An object of class "family" for generalized student distribution.
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
#' elliptical.fitGs <- elliptical(y ~ x1+x2+x3, family = Gstudent(c(1,1))
#' ,data=luz)
#' family(elliptical.fitGs)
#' @export


Gstudent <- function (parm = stop("no s or r argument")) 
{
  if ((parm[1] <= 0) || (parm[2] <= 0)) 
    stop(paste("s and r must be positive"))
  make.family.elliptical("Gstudent", arg = list(s = parm[1], 
                                                r = parm[2]))
}

