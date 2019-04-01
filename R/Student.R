#' @title t-Student Distribution
#' @description Family objects for t-Student distribution provide a convenient way to specify the details of the models used by functions such as \code{elliptical}.
#' @param df degrees of freedom.
#' @return An object of class "family" for student distribution.
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
#' elliptical.fitt <- elliptical(y ~ x1+x2+x3, family = Student(df=5)
#' ,data=luz)
#' family(elliptical.fitt)
#' @export

Student <- function (df = stop("no df argument")) 
{
  if (df < 0) 
    stop(paste("allowed values for degrees of freedom positive"))
  make.family.elliptical("Student", arg = df)
}
