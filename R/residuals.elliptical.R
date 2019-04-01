#' @title Extract Model Residuals
#' @method residuals elliptical
#' @description residuals is a generic function which extracts model residuals from objects returned by modeling functions.
#' @param object fit object for elliptical regression model.
#' @param type an character string that indicates the type of residuals.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return Residuals extracted from the \code{object}.
#' @references Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
#' symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
#' \url{https://doi.org/10.1016/j.spl.2005.04.033}
#' @seealso \code{\link{elliptical}}
#' @keywords elliptical
#' @keywords residuals
#' @examples
#' data(luzdat)
#' y <- luzdat$y
#' x1 <- luzdat$x1 ; x1 <- factor(x1) ; x1 <- C(x1,treatment)
#' x2 <- luzdat$x2
#' x3 <- (luzdat$x2)^2
#' luz <- data.frame(y,x1,x2,x3)
#' elliptical.fitt <- elliptical(y ~ x1+x2+x3, family = Student(df=5)
#' ,data=luz)
#' residuals(elliptical.fitt)
#' @export

residuals.elliptical<- function (object, type = c("stand", "pearson", "response"),...) 
{
  type <- match.arg(type)
  rr <- switch(type, pearson = object$resid/sqrt(object$scalevariance), 
               stand = {
                 Xd <- as.matrix(object$Xmodel)
                 Xdi <- solve(t(Xd) %*% Xd)
                 H <- Xd %*% Xdi %*% t(Xd)
                 H1 <- (1/(object$scalevariance * object$scale)) * 
                   H
                 varr <- object$scalevariance * object$dispersion * 
                   (1 - diag(H1))
                 ri <- object$y - object$fitted
                 ri/sqrt(varr)
               }, response = object$y - object$fitted)
  if (is.null(object$na.action)) 
    rr
  else naresid(object$na.action, rr)
}
