#' @title Extract Residuals for Elliptical Regression Model Fits
#' @method residuals elliptical
#' @description This function compute different types of residuals to the fitted elliptical regression model.
#' @param object an object with the result of the fitted elliptical regression model.
#' @param type a character string that indicates the type of residuals. If is \code{stand} will be computed the standard residuals. 
#' If is \code{ordinal} will be computed the ordinal residuals. If is \code{response} will be computed the response residuals. 
#' If is \code{pearson} will be computed the Pearson residuals. If is \code{desvio} will be computed the desviance residuals.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return Residuals of the specific \code{type} extracted from the \code{object}.
#' @references Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
#' symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
#' \doi{10.1016/j.spl.2005.04.033}
#' @seealso \code{\link{residuals}}, \code{\link{elliptical}}, \code{\link{family.elliptical}}
#' @keywords Elliptical regression models
#' @keywords Residuals
#' @examples
#' data(luzdat)
#' y <- luzdat$y
#' x1 <- luzdat$x1 ; x1 <- factor(x1) ; x1 <- C(x1,treatment)
#' x2 <- luzdat$x2
#' x3 <- (luzdat$x2)^2
#' luz <- data.frame(y,x1,x2,x3)
#' elliptical.fitt <- elliptical(y ~ x1+x2+x3, family = Student(df=5)
#' ,data=luz)
#' residuals(elliptical.fitt, type = "stand")
#' @export

residuals.elliptical<- function (object, type = c("stand", "ordinal", "response", "pearson", "desvio"), ...) 
{
  type <- match.arg(type)

  family <- object$family
  dispersion <-object$dispersion
  rord <- object$y - object$fitted
  resid <- (object$y - object$fitted)/sqrt(dispersion)
    
  Xd <- as.matrix(object$Xmodel)
  Xdi <- solve(t(Xd) %*% Xd)
  H <- Xd %*% Xdi %*% t(Xd)
  H1 <- (1/(object$scalevariance * object$scale)) * H
  varr <- object$scalevariance * (1 - diag(H1))
  rstand <- resid/sqrt(varr)
  
  rpearson <- resid/sqrt(object$scalevariance)
  
  n <- length(resid) ;rdesvio <- rep(0,n)
  u <- (object$y - object$fitted)^2/dispersion
  for(i in 1:n){
    logg0 <- family$g0(0, df = family$df, s = family$s, 
                       r = family$r, alpha = family$alpha, mp = family$mp, 
                       epsi = family$epsi, sigmap = family$sigmap, k = family$k)
    loggu <- family$g0(u[i], df = family$df, s = family$s, 
                       r = family$r, alpha = family$alpha, mp = family$mp, 
                       epsi = family$epsi, sigmap = family$sigmap, k = family$k)
    rdesvio[i] <- (sign(resid[i]))*(2*logg0 - 2*loggu)^(.5)
  }
  
  rr <- switch(type, ordinal = rord, response = resid, pearson = rpearson, stand = rstand, desvio = rdesvio)
  if (is.null(object$na.action)) 
    rr
  else naresid(object$na.action, rr)
}
