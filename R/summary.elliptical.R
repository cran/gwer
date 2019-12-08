#' @title Summarizing Elliptical Model Fits.
#' @method summary elliptical
#' @name summary.elliptical
#' @aliases print.summary.elliptical
#' @description This function produce result summary of the result of the fitted elliptical regression model.
#' @param object an object with the result of the fitted elliptical regression model.
#' @param correlation a logical value to return the correlation for the estimated parameters. If \code{FALSE} (by default) not return the correlation matrix.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return returns an object of class \dQuote{summary.elliptical}, a list with follow components: 
#' \item{coefficients}{the matrix of coefficients, standard errors and significance values for parameters hypothesis test.}
#' \item{dispersion}{either the supplied argument or the estimated dispersion with standard error.}
#' \item{residuals}{the residuals from \code{object}.}
#' \item{cov.unscaled}{the unscaled (dispersion = 1) estimated covariance matrix of the estimated coefficients.}
#' \item{corrrelation}{the matrix of correlation for the estimated parameters.}  
#' \item{family}{family from \code{object}.}
#' \item{loglik}{the likelihood logarithm value from \code{object}.}
#' \item{terms}{the \code{terms} object used.}
#' \item{df}{degrees of fredom from \code{object}.}
#' \item{inter}{number of iterations of optimization process.}
#' \item{nas}{a logical vector indicating if there is \code{na} in estimation of coefficients.}
#' \item{call}{the matched call from \code{object}.}
#' \item{scale}{values of the \code{4d_g} for the specified distribution from \code{object}.}
#' \item{scaledispersion}{values of the \code{4f_g} for the specified distribution from \code{object}.}
#' @references Cysneiros, F. J. A., Paula, G. A., and Galea, M. (2007). Heteroscedastic 
#' symmetrical linear models. Statistics & probability letters, 77(11), 1084-1090. 
#' \url{https://doi.org/10.1016/j.spl.2007.01.012} 
#' @seealso \code{\link{glm}}, \code{\link{elliptical}}, \code{\link{family.elliptical}}
#' @keywords Elliptical models
#' @keywords Linear regression models
#' @examples
#' data(luzdat)
#' y <- luzdat$y
#' x1 <- luzdat$x1 ; x1 <- factor(x1) ; x1 <- C(x1,treatment)
#' x2 <- luzdat$x2
#' x3 <- (luzdat$x2)^2
#' luz <- data.frame(y,x1,x2,x3)
#' elliptical.fitt <- elliptical(y ~ x1+x2+x3, family = Student(df=5)
#' ,data=luz)
#' summary(elliptical.fitt)
#' @rdname summary.elliptical
#' @export


summary.elliptical<-function(object, correlation = TRUE,...) 
{
  coef <- object$coef
  disp <- object$dispersion
  scale <- object$scale
  scaledispersion <- object$scaledispersion
  fixed <- object$fixed
  resid <- object$residuals
  wt <- object$weights
  nas <- is.na(coef)
  n <- length(resid) - sum(wt == 0)
  p <- object$rank
  if (is.null(p)) 
    p <- sum(!nas)
  if (!p) {
    warning("\n This model has zero rank --- no summary is provided")
    return(object)
  }
  rdf <- object$df.resid
  if (is.null(rdf)) 
    rdf <- n - p - sum(wt == 0)
  
  R <- object$R[(1:p), (1:p)]
  Rnames <- dimnames(R)
  covun <- solve(qr(R))

  dimnames(covun) <- Rnames
  rowlen <- sqrt(diag(covun))
  cnames <- names(coef[!nas])
  coef <- matrix(rep(coef[!nas], 4), ncol = 4)
  dimnames(coef) <- list(cnames, c("Value", "Std. Error", 
                                   "z-value", "p-value"))
  coef[, 2] <- rowlen[1:p] %o% sqrt(disp/scale)
  coef[, 3] <- coef[, 1]/coef[, 2]
  coef[, 4] <- 2 * pnorm(-abs(coef[, 3]))
  if (!fixed) {
    disp <- matrix(c(disp, sqrt((4 * disp^2)/(n * scaledispersion))), 
                   ncol = 2)
    dimnames(disp) <- list("dispersion", c("Value", "Std. Error"))
  }
  if (correlation) {
    correl <- covun * outer(1/rowlen, 1/rowlen)
    dimnames(correl) <- Rnames
  }
  else correl <- NULL

  summary <- list(coefficients = coef, dispersion = disp, 
                  fixed = fixed, residuals = resid, cov.unscaled = covun[(1:p),(1:p)], 
                  correlation = correl[(1:p), (1:p)], family = object$family, 
                  loglik = object$loglik, terms = object$terms, df = c(p, rdf, n), 
                  iter = object$iter, nas = nas, call = object$call, scale = scale,
                  scaledispersion = scaledispersion)
#  class(summary) <- "summary"
  attr(summary, "class") <- c("summary.elliptical")
  summary
}

#' @rdname summary.elliptical
#' @method print summary.elliptical
#' @noRd
#' @export
print.summary.elliptical <- function (x, digits = 6, quote = T, ...) 
{
  nas <- x$nas
  p <- sum(!nas)
  coef <- x$coef
  correl <- x$correl
  if (any(nas)) {
    nc <- length(nas)
    cnames <- names(nas)
    coef1 <- array(NA, c(nc, 3), list(cnames, dimnames(coef)[[2]]))
    coef1[!nas, ] <- coef
    coef <- coef1
    if (!is.null(correl)) {
      correl1 <- matrix(NA, nc, nc, dimnames = list(cnames, 
                                                    cnames))
      correl1[!nas, !nas] <- correl[1:p, 1:p]
      correl <- correl1
    }
  }
  if (is.null(digits)) 
    digits <- options()$digits
  else {
    old.digits <- options(digits = digits)
    on.exit(options(old.digits))
  }
  cat("Call: ")
  dput(x$call)
  if (any(nas)) 
    cat("\nCoefficients: (", sum(nas), " not defined because of singularities)\n", 
        sep = "")
  else cat("\nCoefficients:\n")
  print(coef, digits = digits)
  cat(paste("\nScale parameter for", x$family$family, ": "))
  cat(signif(x$dispersion[1], digits = digits), " (", if (x$fixed) 
    "fixed"
    else signif(x$dispersion[2], digits = digits), ")\n")
  int <- attr(x$terms, "intercept")
  if (is.null(int)) 
    int <- 1
  df <- x$df
  nobs <- df[3]
  cat("\nDegrees of Freedom:", nobs, "Total;", x$df[2], "Residual\n")
  cat("-2*Log-Likelihood", format(-2 * x$loglik), "\n")
  cat("\nNumber  Iterations:", format(trunc(x$iter)), "\n")
  if (!is.null(correl)) {
    p <- dim(correl)[2]
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      ll <- lower.tri(correl)
      correl[ll] <- format(round(correl[ll], digits))
      correl[!ll] <- ""
      print(correl[-1, -p, drop = F], quote = F, digits = digits)
    }
  }
  invisible(x)
}



