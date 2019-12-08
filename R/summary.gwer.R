#' @title Summarizing Geographically Weighted Elliptical Regression Model Fits.
#' @method summary gwer
#' @name summary.gwer
#' @aliases summary.gwer.print
#' @description This function produce result summary of the result of the fitted geographically weighted elliptical regression model.
#' @param object an object with the result of the fitted geographically weighted elliptical regression model.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return returns an object of class \dQuote{summary.gwer}, a list with follow components: 
#' \item{coefficients}{the matrix of summarizing coefficients, standard errors and significance values for parameters hypothesis test.}
#' \item{dispersion}{either the supplied argument or the estimated dispersion with standard error.}
#' \item{residuals}{the residuals from \code{object}.}
#' \item{family}{family from \code{object}.}
#' \item{results}{a list of results values for fitted geographically weighted elliptical model.}  
#' \item{spdisp}{a logical value indicating whether the dispersion varies geographically from \code{object}.}
#' \item{df}{degrees of fredom from \code{object}.}
#' \item{terms}{the \code{terms} object used.}
#' \item{inter}{number of iterations of optimization process.}
#' \item{nas}{a logical vector indicating if there is \code{na} in estimation of coefficients.}
#' \item{type}{a character string indicating the type of residuals was obtained from \code{object}}
#' \item{hatmatrix}{a logical value indicating if hatmatrix was obtained from \code{object}}
#' \item{call}{the matched call from \code{object}.}
#' \item{scale}{values of the \code{4d_g} for the specified distribution from \code{object}.}
#' \item{scaledispersion}{values of the \code{4f_g} for the specified distribution from \code{object}.}
#' \item{scalevariance}{values of the scale variance for the specified distribution from \code{object}.}
#' @references Brunsdon, C., Fotheringham, A. S. and Charlton, M. E. (1996). 
#' Geographically weighted regression: a method for exploring spatial nonstationarity.
#' Geographical analysis, 28(4), 281-298. \url{https://doi.org/10.1111/j.1538-4632.1996.tb00936.x}
#' @references Cysneiros, F. J. A., Paula, G. A., and Galea, M. (2007). Heteroscedastic 
#' symmetrical linear models. Statistics & probability letters, 77(11), 1084-1090. 
#' \url{https://doi.org/10.1016/j.spl.2007.01.012} 
#' @seealso \code{\link{gwer}}, \code{\link{gwer.sel}}, \code{\link{family.elliptical}}
#' @keywords Geographically weighted regression
#' @keywords Elliptical model
#' @examples
#' data(columbus, package="spData")
#' fit.lm <- lm(CRIME ~ INC, data=columbus)
#' summary(fit.lm)
#' gwer.bw <- gwer.sel(CRIME ~ INC, data=columbus, family = Normal(),
#'                  coords=cbind(columbus$X, columbus$Y))
#' fit.gwer <- gwer(CRIME ~ INC, family = Normal(), bandwidth = gwer.bw, hatmatrix = TRUE,
#'                  spdisp = TRUE, parplot = TRUE, data=columbus, method = "gwer.fit",
#'                  coords=cbind(columbus$X, columbus$Y))
#' summary(fit.gwer)  
#' \donttest{
#' data(columbus, package="spData")
#' fit.elliptical <- elliptical(CRIME ~ INC, family = Student(df=4), data=columbus)
#' summary(fit.elliptical)
#' gwer.bw <- gwer.sel(CRIME ~ INC, data=columbus, family = Student(df=4),
#'                  coords=cbind(columbus$X, columbus$Y), method = 'aic')
#' gwer.fitt <- gwer(CRIME ~ INC, family = Student(df=4), bandwidth = gwer.bw, hatmatrix = TRUE,
#'                  spdisp = TRUE, parplot = TRUE, data=columbus, method = "gwer.fit",
#'                  coords=cbind(columbus$X, columbus$Y))
#' summary(gwer.fitt)  
#' }
#' @rdname summary.gwer
#' @export


summary.gwer<-function(object, ...) 
{
  family <- object$family
  coef <- object$coef
  disp <- object$dispersion
  scale <- 4 * family$g2(resid, df = family$df, 
                         r = family$r, s = family$s, alpha = family$alpha, 
                         mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                         k = family$k)    
  scaledispersion <- -1 + 4 * family$g3(resid, 
                                        df = family$df, r = family$r, s = family$s, alpha = family$alpha, 
                                        mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                                        k = family$k)
  scalevariance <- family$g4(resid, df = family$df, 
                             r = family$r, s = family$s, alpha = family$alpha, 
                             mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                             k = family$k)
  fixed <- object$lm$fixed
  resid <- residuals.gwer(object, type = object$type)
  wt <- object$lm$weights
  nas <- is.na(coef)
  n <- length(resid) - sum(wt == 0)
  p <- object$lm$rank
  if (is.null(p)) 
    p <- sum(!nas)
  if (!p) {
    warning("\n This model has zero rank --- no summary is provided")
    return(object)
  }
  rdf <- object$lm$df.resid
  if (is.null(rdf)) 
    rdf <- n - p - sum(wt == 0)

#  R <- object$R[(1:p), (1:p)]
#  Rnames <- dimnames(R)
#  covun <- solve(qr(R))
#  dimnames(covun) <- Rnames
#  rowlen <- sqrt(diag(covun))
#  cnames <- names(coef[!nas])
  resid<- summary(resid) ; CM <- list() 
  coefest <- as.data.frame(object$coef$est)[, , drop=FALSE]
  if (any(is.na(coefest))) {
    coefest <- na.omit(coefest)
    warning("NAs in coefficients dropped")
  }
  CM$est <- t(apply(coefest, 2, summary))
  if (is.null(dim(CM$est))) CM$est <- t(as.matrix(CM$est))
  if (!object$fp.given) {
    if(object$spdisp){
      CM$est <- cbind(CM$est, as.vector(c(object$lm$coefficients,object$lm$dispersion)))
    } else {
      CM$est <- cbind(CM$est, as.vector(object$lm$coefficients))
    }
    colnames(CM$est) <- c(colnames(CM$est)[1:6], "Global")
  }
  
  if(!object$fp.given && object$hatmatrix){
    coefse <- as.data.frame(object$coef$se)[, , drop=FALSE]
    coefpvalue <- as.data.frame(object$coef$pvalue)[, , drop=FALSE]
    if (any(is.na(coefse)) || any(is.na(coefpvalue))) {
      coefse <- na.omit(coefse)
      coefpvalue <- na.omit(coefpvalue)
      warning("NAs in statistcs dropped")
    }
    CM$se <- t(apply(coefse, 2, summary))
    if (is.null(dim(CM$se))) CM$se <- t(as.matrix(CM$se))
    if (!object$fp.given) {
      if(object$spdisp){
        CM$se <- cbind(CM$se, as.vector(c(summary(object$lm)$coefficients[,2],summary(object$lm)$dispersion[2])))
      } else {
        CM$se <- cbind(CM$se, as.vector(summary(object$lm)$coefficients[,2]))
      }
      colnames(CM$se) <- c(colnames(CM$se)[1:6], "Global")
    }
    CM$pvalue <- t(apply(coefpvalue, 2, summary))
    if (is.null(dim(CM$pvalue))) CM$pvalue <- t(as.matrix(CM$pvalue))
    if (!object$fp.given) {
      CM$pvalue <- cbind(CM$pvalue, as.vector(summary(object$lm)$coefficients[,4]))
      colnames(CM$pvalue) <- c(colnames(CM$pvalue)[1:6], "Global")
    }
  }
#  dimnames(coef) <- list(cnames, c("Value", "Std. Error", 
#                                   "z-value", "p-value"))
#  coef[, 2] <- rowlen[1:p] %o% sqrt(disp/scale)
#  coef[, 3] <- coef[, 1]/coef[, 2]
#  coef[, 4] <- 2 * pnorm(-abs(coef[, 3]))
  if (!fixed && !object$spdisp) {
    disp <- matrix(c(disp, sqrt((4 * disp^2)/(n * scaledispersion))), 
                   ncol = 2)
    dimnames(disp) <- list("dispersion", c("Value", "Std. Error"))
  }
#  if (correlation) {
#    correl <- covun * outer(1/rowlen, 1/rowlen)
#    dimnames(correl) <- Rnames
#  }
#  else correl <- NULL

  summary <- list(coefficients = CM, dispersion = disp, fixed = fixed, residuals = resid, 
                  #cov.unscaled = covun[(1:p),(1:p)], correlation = correl[(1:p), (1:p)], 
                  family = object$family, results = object$results, spdisp = object$spdisp, 
                  df = c(p, rdf, n), terms = object$terms, nas = nas, type = object$type,
                  hatmatrix = object$hatmatrix, call = object$this.call, scale = scale, 
                  scaledispersion = scaledispersion, scalevariance = scalevariance)

  attr(summary, "class") <- c("summary.gwer")
  summary
}

#' @rdname summary.gwer
#' @method print summary.gwer
#' @noRd
#' @export
print.summary.gwer <- function (x, digits = 6, ...) 
{
  nas <- x$nas
  p <- sum(!nas)
  coef <- x$coefficients
#  correl <- x$correl
  if (any(nas)) {
    nc <- length(nas)
    cnames <- names(nas)
    coef1 <- array(NA, c(nc, 3), list(cnames, dimnames(coef)[[2]]))
    coef1[!nas, ] <- coef
    coef <- coef1
#    if (!is.null(correl)) {
#      correl1 <- matrix(NA, nc, nc, dimnames = list(cnames, 
#                                                    cnames))
#      correl1[!nas, !nas] <- correl[1:p, 1:p]
#      correl <- correl1
#    }
  }
  if (is.null(digits)) 
    digits <- options()$digits
  else {
    old.digits <- options(digits = digits)
    on.exit(options(old.digits))
  }
  cat("Call: ")
  dput(x$call)
  
  cat("\nSummary of ", x$type ," Residuals:\n")
  print(x$residuals, digits = digits)
  if (any(nas)) 
    cat("\nCoefficients: (", sum(nas), " not defined because of singularities)\n", 
        sep = "")
  else cat("\nSummary of GWER coefficient estimates:\n")
  print(coef$est, digits = digits)
  if(!is.null(coef$se) && !is.null(coef$pvalue)){
    cat("\nSummary of GWER standard errors estimates:\n")
    print(coef$se, digits = digits)
#    cat("\nSummary of GWER p-Values for significance test:\n")
#    print(coef$pvalue, digits = digits)   
  }
  
  cat(paste("\n Global scale parameter for", x$family$family, ": "))
  cat(signif(x$dispersion[1], digits = digits), " (", if (x$fixed) 
    "fixed"
    else signif(x$dispersion[2], digits = digits), ")\n")
  int <- attr(x$terms, "intercept")
  if (is.null(int)) 
    int <- 1
  df <- x$df ;  nobs <- df[3]
#  cat("\nDegrees of Freedom:", nobs, "Total;", x$df[2], "Residual\n")
  cat("Number of data points:", nobs, "\n")
#  cat("\nNumber  Iterations:", format(trunc(x$iter)), "\n")
#  if(!x$spdisp)
#    cat("Global dispersion:", x$dispersion, "\n")
  if (x$hatmatrix) {
    cat("Effective number of parameters (residual: 2traceS - traceS'S):", 2*x$results$nu1 -
          x$results$nu2, "\n")
    cat("Effective degrees of freedom (residual: 2traceS - traceS'S):", x$results$edf, "\n")
    cat("Effective number of parameters (model: traceS):",
        x$results$nu1, "\n")
    cat("Effective degrees of freedom (model: traceS):",
        (nobs - x$results$nu1), "\n")
#    cat("-2*Log-Likelihood:", format(x$results$logLik), "\n")
    cat("-2*Log-Likelihood:", format(x$results$logLik), " AIC:", x$results$AIC, " AICc:", 
        x$results$AICc, " BIC:", x$results$BIC, "\n")
  }
  #if (!is.null(correl)) {
  #  p <- dim(correl)[2]
  #  if (p > 1) {
  #    cat("\nCorrelation of Coefficients:\n")
  #    ll <- lower.tri(correl)
  #    correl[ll] <- format(round(correl[ll], digits))
  #    correl[!ll] <- ""
  #    print(correl[-1, -p, drop = F], quote = F, digits = digits)
  #  }
  #}
  invisible(x)
}



