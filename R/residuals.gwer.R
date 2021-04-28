#' @title Extract Residuals for Geographically Weighted Elliptical Regression Model Fits
#' @method residuals gwer
#' @description This function compute different types of residuals to the fitted geographically weighted elliptical regression model.
#' @param object an object with the result of the fitted geographically weighted elliptical regression model.
#' @param type a character string that indicates the type of residuals. If is \code{stand} will be computed the standar residuals. 
#' If is \code{ordinal} will be computed the ordinal residuals. If is \code{response} will be computed the response residuals. 
#' If is \code{pearson} will be computed the pearson residuals. If is \code{desvio} will be computed the desviance residuals.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return Residuals of the specific \code{type} extracted from the \code{object}.
#' @references Brunsdon, C., Fotheringham, A. S. and Charlton, M. E. (1996). 
#' Geographically weighted regression: a method for exploring spatial nonstationarity.
#' Geographical analysis, 28(4), 281-298. \doi{10.1111/j.1538-4632.1996.tb00936.x}
#' @references Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
#' symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
#' \doi{10.1016/j.spl.2005.04.033}
#' @seealso \code{\link{residuals}}, \code{\link{gwer}}, \code{\link{family.elliptical}}
#' @keywords Geographically weighted regression
#' @keywords Elliptical regression models
#' @keywords Residuals
#' @examples
#' \donttest{
#' data(georgia, package = "spgwr")
#' fit.formula <- PctBach ~ TotPop90 + PctRural + PctFB + PctPov
#' gwer.bw.t <- bw.gwer(fit.formula, data = gSRDF, family = Student(3), adapt = TRUE)
#' gwer.fit.t <- gwer(fit.formula, data = gSRDF, family = Student(3), bandwidth = gwer.bw.t, 
#'                    adapt = TRUE, parplot = FALSE, hatmatrix = TRUE, spdisp = TRUE, 
#'                    method = "gwer.fit")
#' summary(gwer.fit.t) 
#' residuals(gwer.fit.t, type = "stand") 
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


#  if (charmatch(dist, "Normal", F)) {
#    dist.y <- pnorm(object$y, object$fitted, dispersion)
#  }
#  else if (charmatch(dist, "Cauchy", F)) {
#    dist.y <- pcauchy(object$y, object$fitted, dispersion)
#  }
#  else if (charmatch(dist, "Student", F)) {
#    dist.y <- pt(resid, family$df)
#  }
#  else if (charmatch(dist, "Gstudent", F)) {
#    dist.y <- pgstudent(resid, family$s, family$r)
#  }
#  else if (charmatch(dist, "LogisI", F)) {
#    stop(paste("not implemented yet"))
#    dist.y <- plogisI(object$y, object$fitted, dispersion)
#  }
#  else if (charmatch(dist, "LogisII", F)) {
#    dist.y <- plogisII(resid)
#  }
#  else if (charmatch(dist, "Glogis", F)) {
#    stop(paste("not implement yet"))
#    dist.y <- pglogis(resid, family$alpha, family$mp)
#  }
#  else if (charmatch(dist, "Cnormal", F)) {
#    stop(paste("not implemented yet"))
#    dist.y <- NULL#pcnormal(resid, family$epsi, family$sigmap)
#  }
#  else if (charmatch(dist, "Powerexp", F)) {
#    dist.y <- ppowerexp(resid, family$k)
#  }
#  rquant <- qnorm(dist.y, 0, 1)
  
  rr <- switch(type, ordinal = rord, response = resid, pearson = rpearson, stand = rstand, desvio = rdesvio)#, quantile = rquant)
  if (is.null(object$na.action)) 
    rr
  else naresid(object$na.action, rr)
}
