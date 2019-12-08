#' @title Extract Residuals for Geographically Weighted Elliptical Regression Model Fits
#' @method residuals gwer
#' @description This function compute differents type of residuals to the fittedgeographically weighted elliptical regression model.
#' @param object an object with the result of the fitted geographically weighted elliptical regression model.
#' @param type a character string that indicates the type of residuals. If is \code{stand} will be computed the standar residuals. 
#' If is \code{ordinal} will be computed the ordinal residuals. If is \code{response} will be computed the response residuals. 
#' If is \code{pearson} will be computed the pearson residuals. If is \code{desvio} will be computed the desviance residuals.
#' By default is \code{stand}.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return Residuals of the specific \code{type} extracted from the \code{object}.
#' @references Brunsdon, C., Fotheringham, A. S. and Charlton, M. E. (1996). 
#' Geographically weighted regression: a method for exploring spatial nonstationarity.
#' Geographical analysis, 28(4), 281-298. \url{https://doi.org/10.1111/j.1538-4632.1996.tb00936.x}
#' @references Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
#' symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
#' \url{https://doi.org/10.1016/j.spl.2005.04.033}
#' @seealso \code{\link{residuals}}, \code{\link{gwer}}, \code{\link{family.elliptical}}
#' @keywords Geographically weighted regression
#' @keywords Elliptical models
#' @keywords Residuals
#' @examples
#' data(columbus, package="spData")
#' fit.lm <- lm(CRIME ~ INC, data=columbus)
#' summary(fit.lm)
#' gwer.bw <- gwer.sel(CRIME ~ INC, data=columbus, family = Normal(),
#'                  coords=cbind(columbus$X, columbus$Y))
#' fit.gwer <- gwer(CRIME ~ INC, family = Normal(), bandwidth = gwer.bw, 
#'                  spdisp = TRUE, parplot = TRUE, data=columbus, method = "gwer.fit",
#'                  coords=cbind(columbus$X, columbus$Y))
#' residuals(fit.gwer, type = "stand")
#' \donttest{
#' data(columbus, package="spData")
#' fit.elliptical <- elliptical(CRIME ~ INC, family = Student(df=4), data=columbus)
#' summary(fit.elliptical)
#' gwer.bw <- gwer.sel(CRIME ~ INC, data=columbus, family = Student(df=4),
#'                  coords=cbind(columbus$X, columbus$Y), method = 'aic')
#' gwer.fitt <- gwer(CRIME ~ INC, family = Student(df=4), bandwidth = gwer.bw, hatmatrix = TRUE,
#'                  spdisp = TRUE, parplot = TRUE, data=columbus, method = "gwer.fit",
#'                  coords=cbind(columbus$X, columbus$Y))
#' residuals(gwer.fitt)  
#' }
#' @export

residuals.gwer <- function (object, type = c("stand", "ordinal", "response", "pearson", "desvio"), ...) 
{
  type <- match.arg(type)

  family <- object$family
  dispersion <- object$dispersion
  if(class(object)[1]=="gwer"){
    object$y <- object$lm$y
    object$Xmodel <- object$lm$Xmodel
  }  
  rord <- (object$y - object$fitted)
  resid <- (object$y - object$fitted)/sqrt(dispersion)
  n <- length(resid) ;rdesvio <- rep(0,n)
  u <- (object$y - object$fitted)^2/dispersion
  
  Xd <- as.matrix(object$Xmodel)
  H <- matrix(0,n,n)
  if(class(object)[1]=="gwer"){
    Wi <- object$gweights
    for(i in 1:n){
      H[i,] <- Xd[i,] %*% solve(t(Xd) %*% diag(Wi[i,]) %*% Xd) %*% t(Xd) %*% diag(Wi[i,])  
    }
  } else {
    Wi <- diag(object$gweights)
    H <- Xd %*% solve(t(Xd) %*% Wi %*% Xd) %*% t(Xd) %*% Wi
  }

    
  scale <- 4 * family$g2(resid, df = family$df, 
                                r = family$r, s = family$s, alpha = family$alpha, 
                                mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                                k = family$k) 
  scalevariance <- family$g4(resid, df = family$df, 
                             r = family$r, s = family$s, alpha = family$alpha, 
                             mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                             k = family$k)
  
  H1 <- (1/(scalevariance * scale)) * H
  varr <- scalevariance * (1 - diag(H1))
  rstand <- resid/sqrt(varr)
  
  rpearson <- resid/sqrt(scalevariance)
  

  for(i in 1:n){
    logg0 <- family$g0(0, df = family$df, s = family$s, 
                       r = family$r, alpha = family$alpha, mp = family$mp, 
                       epsi = family$epsi, sigmap = family$sigmap, k = family$k)
    loggu <- family$g0(u[i], df = family$df, s = family$s, 
                       r = family$r, alpha = family$alpha, mp = family$mp, 
                       epsi = family$epsi, sigmap = family$sigmap, k = family$k)
    rdesvio[i] <- sqrt(Wi[i,i])*(sign(resid[i]))*(2*logg0 - 2*loggu)^(.5)
  }

  rres <- object$y - object$fitted
  
  rr <- switch(type, ordinal = rord, response = resid, pearson = rpearson, stand = rstand, desvio = rdesvio)
  if (is.null(object$na.action)) 
    rr
  else naresid(object$na.action, rr)
}
