#' @title Simulated Envelope of Residuals for Geographically Weighted Elliptical Regression Models
#' @import glogis
#' @description This function produces quantile-quantile residuals plot with simulated envelope for the specified error 
#' distribution in geographically weighted elliptical regression Models.
#' @param object an object with the result of the fitted geographically weighted elliptical regression Models.
#' @param B number of Monte Carlo simulations.
#' @param arg a numerical or vector representing the distribution parameters used in fitted model.
#' @param xlab a label for the x axis, defaults to a description of x.
#' @param ylab a label for the y axis, defaults to a description of y.
#' @param ident a numerical indicate the number of observation identified in plot.
#' @param ident.labels an optional character vector giving labels for the identified points.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @references Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
#' symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
#' \doi{10.1016/j.spl.2005.04.033}
#' @seealso \code{\link{glm}}, \code{\link{elliptical}}, \code{\link{family.elliptical}}
#' @keywords Geographically Weighted Regression
#' @keywords Elliptical regression models
#' @keywords quantile-quantile plots
#' @examples
#' \donttest{
#' data(georgia, package = "spgwr")
#' fit.formula <- PctBach ~ TotPop90 + PctRural + PctFB + PctPov
#' gwer.bw.t <- bw.gwer(fit.formula, data = gSRDF, family = Student(3), adapt = TRUE)
#' gwer.fit.t <- gwer(fit.formula, data = gSRDF, family = Student(3), bandwidth = gwer.bw.t, 
#'                    adapt = TRUE, parplot = FALSE, hatmatrix = TRUE, spdisp = TRUE, 
#'                    method = "gwer.fit")
#' }
#' @export

gwer.envelope <- function (object, B = 100, arg, xlab = NULL, ylab = NULL, ident = NULL, ident.labels = NULL, ...) 
{
  initial <- NULL
  #data <- object$lm$data
  if(is.null(xlab))
    xlab <- "Quantiles of N(0,1)"
  if(is.null(ylab))
    ylab <- "Standardized residual"
  X <- model.matrix(object$lm$terms)
  Xd <- as.matrix(object$lm$Xmodel)
  n <- nrow(Xd)
  p <- ncol(Xd)
  W <- object$gweights
  l.fitted = object$flm
  family <- object$family
  control <- object$lm$control
#  ro <- object$resid
#  tdf <- ro/sqrt(object$scalevariance)
  resid <- (object$y - object$fitted)/sqrt(object$dispersion)
  H <- matrix(0,n,n)
  for(i in 1:n)
      H[i,] <- Xd[i,] %*% solve(t(Xd) %*% diag(W[i,]) %*% Xd) %*% t(Xd) %*% diag(W[i,])  
  
  scale <- 4 * family$g2(resid, df = family$df, 
                         r = family$r, s = family$s, alpha = family$alpha, 
                         mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                         k = family$k) 
  scalevariance <- family$g4(resid, df = family$df, 
                             r = family$r, s = family$s, alpha = family$alpha, 
                             mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                             k = family$k)
  H1 <- (1/(scalevariance * scale)) * H
  varr <- scalevariance * object$dispersion * (diag(1, n) - H1)
  varr <- diag(varr)
  ri <- object$lm$y - object$fitted
  tdf <- ri/sqrt(varr)
  e <- e.i <- matrix(0, n, B) ; med <- matrix(0, n, n)
  resp <- NULL
  method <- "gwer.fit"
  elliptical.fitter <- get(method)
#  offset = object$offset
#  if (length(offset) == 1 && offset == 0) 
#    offset <- rep(0, nobs)

  for(i in 1:n){
    mu <- object$flm[i,]
    phi <- object$dispersion[i]
    w.i <- W[i,]
    for (j in 1:B) {
      dist <- object$family[[1]]
      if (charmatch(dist, "Normal", F)) {
        resp <- rnorm(n, 0, 1)
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = Normal(), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "Cauchy", F)) {
        resp <- rcauchy(n, 0, 1)
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = Cauchy(), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "Student", F)) {
        resp <- rt(n, arg)
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = Student(arg), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "Gstudent", F)) {
        resp <- rgstudent(n, arg[1], arg[2])
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = Gstudent(arg), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "LogisI", F)) {
        stop(paste("not implemented yet"))
        resp <- rlogisI(n, 0, 1)
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = LogisI(), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "LogisII", F)) {
        resp <- rlogisII(n)
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = LogisII(), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "Glogis", F)) {
        stop(paste("not implement yet"))
        resp <- rglogis(n, arg[1], arg[2])
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = Glogis(arg), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "Cnormal", F)) {
        stop(paste("not implemented yet"))
        resp <- rcnormal(n, arg[1], arg[2])
        #fit <- elliptical(resp ~ X + (-1), family = Cnormal(arg), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "Powerexp", F)) {
        resp <- rpowerexp(n, arg)
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = Powerexp(arg), 
        #                  control = glm.control(maxit = 1000))
      }

      lm.i <- elliptical.fitter(X = X, Y = resp, gweights = w.i, family = family, offset = NULL,
                                dispersion = NULL, maxit = control$maxit, epsilon = control$epsilon,
                                trace = control$trace, ...)
      ro.i <- lm.i$resid
      td.i <- ro.i/sqrt(lm.i$scalevariance)
      Xd.i <- as.matrix(lm.i$Xmodel)
      H.i <- Xd.i %*% solve(t(Xd.i) %*% diag(w.i) %*% Xd.i) %*% t(Xd.i) %*% diag(w.i) 
      H1.i <- (1/(lm.i$scalevariance * lm.i$scale)) * H.i
      varr.i <- lm.i$scalevariance * lm.i$dispersion * (diag(1, n) - H1.i)
      varr.i <- diag(varr.i)
      ri.i <- resp - lm.i$fitted
      td.i <- ri.i/sqrt(varr.i)
      e.i[, j] <- td.i #sort(td.i)
    }
    e[i, ] <- e.i[i, ]
  }
  e <- apply(e, 2, sort)

  e1 <- numeric(n)
  e2 <- numeric(n)
  e3 <- numeric(n)
  e4 <- numeric(n)
  e5 <- numeric(n)
  e6 <- numeric(n)
  e7 <- numeric(n)
  for (i in 1:n) {
    eo <- sort(e[i, ])
    e1[i] <- eo[ceiling(B * 0.05)]
    e2[i] <- eo[ceiling(B * 0.95)]
  }
  e3 <- t(t(apply(e, 2, mean)))
  e4 <- t(t(apply(e, 2, vari)))
  e5 <- t(t(apply(e, 2, skewn)))
  e6 <- t(t(apply(e, 2, kurt)))
  e7 <- cbind(e3, e4, e5, e6)
  desc <- apply(e, 2, mean)
  med <- apply(e, 1, mean)
  faixa <- range(tdf, e1, e2)
  screen(4)
  par(pty = "s")
  points.p <- qqnorm(tdf, xlab = xlab, ylab = ylab, ylim = faixa, pch = 16, main = '')
  par(new = TRUE)
  qqnorm(e1, axes = F, xlab = "", ylab = "", type = "l", ylim = faixa, 
         lty = 1, main = '')
  par(new = TRUE)
  qqnorm(e2, axes = F, xlab = "", ylab = "", type = "l", ylim = faixa, 
         lty = 1, main = '')
  par(new = TRUE)
  qqnorm(med, axes = F, xlab = "", ylab = "", type = "l", 
         ylim = faixa, lty = 2, main = '')
  if(is.null(ident.labels))
    ident.labels <- seq_along(points.p$x)
  if(!is.null(ident))
    identify(points.p$x, points.p$y, labels = ident.labels, n = ident)
  x <- list(mean = desc[1], var = desc[2], skewness = desc[3], 
            kurtosis = desc[4], tdf = tdf)
  invisible(x)
}



vari <- function (x) 
{
  wnas <- x[!is.na(x)]
  var(x, na.rm = TRUE) * (length(wnas) - 1)/length(wnas)
}


skewn <- function (x, na.rm = F, method = "fisher") 
{
  method <- char.expand(method, c("fisher", "moment"), stop("argument 'method' must match either \"fisher\" or \"moment\""))
  if (na.rm) {
    wnas <- x[!is.na(x)]
    if (length(wnas)) 
      x <- wnas
  }
  else if (any(is.na(x[!is.na(x)]))) 
    return(NA)
  n <- length(x)
  if (method == "fisher" && n < 3) 
    return(NA)
  x <- x - mean(x)
  if (method == "moment") 
    (sum(x^3)/n)/(sum(x^2)/n)^1.5
  else ((sqrt(n * (n - 1))/(n - 2)) * (sum(x^3)/n))/((sum(x^2)/n)^1.5)
}

kurt <- function (x, na.rm = F, method = "fisher") 
{
  method <- char.expand(method, c("fisher", "moment"), stop("argument 'method' must match either \"fisher\" or \"moment\""))
  if (na.rm) {
    wnas <- x[!is.na(x)]
    if (length(wnas)) 
      x <- wnas
  }
  else if (any(is.na(x[!is.na(x)]))) 
    return(NA)
  n <- length(x)
  if (method == "fisher" && n < 4) 
    return(NA)
  x <- x - mean(x)
  if (method == "moment") 
    (sum(x^4)/n)/(sum(x^2)/n)^2 - 3
  else ((n + 1) * (n - 1) * ((sum(x^4)/n)/(sum(x^2)/n)^2 - 
                               (3 * (n - 1))/(n + 1)))/((n - 2) * (n - 3))
}

