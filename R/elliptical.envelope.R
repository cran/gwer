#' @title Simulated Envelope of Residuals for Elliptical Regression Models
#' @import glogis
#' @description This function produces quantile-quantile residuals plot with simulated envelope for the specified error distribution in elliptical regression models.
#' @param object an object with the result of the fitted elliptical regression model.
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
#' @keywords Elliptical regression models
#' @keywords Quantile-quantile plots
#' @examples
#' data(luzdat)
#' y <- luzdat$y
#' x1 <- luzdat$x1 ; x1 <- factor(x1) ; x1 <- C(x1,treatment)
#' x2 <- luzdat$x2
#' x3 <- (luzdat$x2)^2
#' luz <- data.frame(y,x1,x2,x3)
#' elliptical.fitt <- elliptical(y ~ x1+x2+x3, family = Student(df=5),
#' data=luz)
#' elliptical.envelope(elliptical.fitt, B=100, arg=5)
#' @export

elliptical.envelope <- function (object, B = 100, arg, xlab = NULL, ylab = NULL, ident = NULL, ident.labels = NULL, ...) 
{
  initial <- NULL
  if(is.null(xlab))
    xlab <- "Quantiles of N(0,1)"
  if(is.null(ylab))
    ylab <- "Standardized residual"
  X <- model.matrix(object$terms)
  Xd <- as.matrix(object$Xmodel)
  n <- nrow(Xd)
  p <- ncol(Xd)
  family <- object$family
  control <- object$control
#  ro <- object$resid
#  tdf <- ro/sqrt(object$scalevariance)
  H <- Xd %*% solve(t(Xd) %*% Xd) %*% t(Xd)
  H1 <- (1/(object$scalevariance * object$scale)) * H
  varr <- object$scalevariance * object$dispersion * (diag(1, n) - H1)
  varr <- diag(varr)
  ri <- object$y - object$fitted
  tdf <- ri/sqrt(varr)
  e <- matrix(0, n, B)
  mu <- object$fitted
  phi <- object$dispersion
  resp <- NULL
  method <- "elliptical.fit"
  elliptical.fitter <- get(method)
  #  offset = object$offset
  #  if (length(offset) == 1 && offset == 0) 
  #    offset <- rep(0, nobs)
  
  for (i in 1:B) {
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

    fit <- elliptical.fitter(X = X, Y = resp, family = family, dispersion = NULL, offset = 
                             NULL, maxit = control$maxit, epsilon = control$epsilon,
                             trace = control$trace, ...)
    ro <- fit$resid
    td <- ro/sqrt(fit$scalevariance)
    Xd <- as.matrix(fit$Xmodel)
    H <- Xd %*% solve(t(Xd) %*% Xd) %*% t(Xd)
    H1 <- (1/(fit$scalevariance * fit$scale)) * H
    varr <- fit$scalevariance * fit$dispersion * (diag(1, n) - H1)
    varr <- diag(varr)
    ri <- resp - fit$fitted
    td <- ri/sqrt(varr)
    e[, i] <- sort(td)
  }

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
  desc <- apply(e7, 2, mean)
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
            kurtosis = desc[4])
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

