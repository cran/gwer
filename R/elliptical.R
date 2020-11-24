#' @title Elliptical Regression Models
#' @import stats
#' @import assertthat
#' @name elliptical
#' @aliases print.elliptical
#' @description The function implements linear elliptical regression models. This models is specified giving a symbolic description of the systematic and stochastic components.
#' @param formula regression model formula as in \code{glm}.
#' @param family a description of the error distribution to be used in the model (see \code{elliptical.family} for details of family functions).
#' @param data an optional data frame, list or environment containing the variables in the model.
#' @param dispersion an optional fixed value for dispersion parameter.
#' @param weights an optional numeric vector of weights to be used in the fitting process.
#' @param subset an optional numeric vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain NAs (see \code{glm}).
#' @param method optimization method used to estimate the parameters. The default method "elliptical.fit" uses Fisher's scoring method. The alternative "model.frame" returns the model frame and does no fitting.
#' @param control a list of parameters for controlling the fitting process. For \code{elliptical} this is passed by \code{glm.control}.
#' @param model a logical value indicating whether model frame should be included as a component of the return.
#' @param x a logical value indicating whether the response vector used in the fitting process should be returned as components of the return.
#' @param y a logical value indicating whether model matrix used in the fitting process should be returned as components of the return.
#' @param contrasts an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.
#' @param offset this can be used to specify an a priori known component to be included in the linear predictor during fitting as in \code{glm}.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return returns an object of class \dQuote{elliptical}, a list with follow components: 
#' \item{coefficients}{coefficients of location parameters.}
#' \item{dispersion}{coefficient of dispersion parameter.}
#' \item{residuals}{standardized residuals.}
#' \item{fitted.values}{the fitted mean values.}
#' \item{loglik}{the likelihood logarithm value for the fitted model.}  
#' \item{Wg}{values of the function \code{W_g(u)}.}
#' \item{Wgder}{values for the function \code{W^{(1)}_g(u)}.}
#' \item{v}{values for the function \code{V(u)}.}
#' \item{rank}{the numeric rank for the fitted model.}
#' \item{R}{the matrix of correlation for the estimated parameters.}
#' \item{inter}{number of iterations of optimization process.}
#' \item{scale}{values of the \code{4d_g} for the specified distribution.}
#' \item{scaledispersion}{values of the \code{4f_g} for the specified distribution.}
#' \item{scalevariance}{values of the scale variance for the specified distribution.}
#' \item{df}{degree of freedom for t-student distribution.}
#' \item{s, r}{shape parameters for generalized t-student distribution.}
#' \item{alpha}{shape parameter for contaminated normal and generalized logistic distributions.}  
#' \item{mp}{shape parameter for generalized logistic distribution.}
#' \item{epsi,sigmap}{dispersion parameters for contaminated normal distribution.}
#' \item{k}{shape parameter for power exponential distribution.}
#' \item{Xmodel}{the model matrix.}
#' \item{weights}{the working weights, that is the weights in the final iteration of optimization process}
#' \item{df.residuals}{the residual degrees of freedom.}
#' \item{family}{the \code{family} object used.}
#' \item{formula}{the \code{formula} supplied.}
#' \item{terms}{the \code{terms} object used.}
#' \item{contrasts}{(where relevant) the contrasts used.}
#' \item{control}{the value of the \code{control} argument used.}
#' \item{call}{the matched call.}
#' \item{y}{the response variable used.}
#' @references Cysneiros, F. J. A., Paula, G. A., and Galea, M. (2007). Heteroscedastic 
#' symmetrical linear models. Statistics & probability letters, 77(11), 1084-1090. 
#' \url{https://doi.org/10.1016/j.spl.2007.01.012} 
#' @references Fang, K. T., Kotz, S. and NG, K. W. (1990, ISBN:9781315897943).
#' Symmetric Multivariate and Related Distributions. London: Chapman and Hall.
#' @seealso \code{\link{glm}}, \code{\link{family.elliptical}}, \code{\link{summary.elliptical}}
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
#' elliptical.fitLII <- elliptical(y ~ x1+x2+x3, family = LogisII()
#' ,data=luz)
#' @rdname elliptical
#' @export

elliptical <- function (formula = formula(data), family = Normal, data, dispersion = NULL, 
                       weights, subset, na.action = "na.fail", method = "elliptical.fit", 
                       control = glm.control(epsilon = 1e-04, maxit = 100, trace = F), 
                       model = F, x = F, y = T, contrasts = NULL, offset, ...) 
{
  call <- match.call()
  dist <- as.character(call$family)[1]
  user.def <- F
  if (charmatch(dist, c("Normal", "Cauchy", "Student", "Gstudent", 
                        "LogisI", "LogisII", "Glogis", "Cnormal", "Powerexp"), 
                nomatch = F)) 
    dist <- match.arg(dist, c("Normal", "Cauchy", "Student", 
                              "Gstudent", "LogisI", "LogisII", "Glogis", "Cnormal", 
                              "Powerexp"))
  else user.def <- T
  if (!charmatch(method, c("model.frame", "elliptical.fit"), 
                 F)) 
    stop(paste("\n unimplemented method:", method))
  if (missing(data)) 
    data <- environment(formula)
  m <- match.call(expand.dots = F)
  m$family <- m$method <- m$control <- m$model <- m$dispersion <- m$x <- m$y <- m$contrasts <-  m$offset <- m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  if (method == "model.frame") 
    return(m)
  {
    if (!missing(family) && !charmatch(dist, c("Normal", 
                                               "Cauchy", "Student", "Gstudent", "LogisI", "LogisII", 
                                               "Glogis", "Cnormal", "Powerexp"), F)) 
      cat(paste("\n work with user-defined family:", call$family, 
                "\n"))
    }
  if (!missing(dispersion) && is.number(dispersion) && !(dispersion > 
                                                         0)) 
    stop("\n no negative values for dispersion parameter")
  Terms <- attr(m, "terms")
  Y <- model.extract(m, "response")
  if (!is.numeric(Y)) 
    stop("\n response must be numeric")
  X <- model.matrix(Terms, m, contrasts)
  if (!is.numeric(X)) 
    stop("\n model matrix must be numeric")
  offset <- model.extract(m, offset)
  nobs <- nrow(X)
  if (length(offset) == 1 && offset == 0) 
    offset <- rep(0, nobs)
  w <- model.extract(m, weights)
  wzero <- rep(F, nrow(m))
  if (!length(w)) 
    w <- rep(1, nrow(m))
  else if (any(w < 0)) 
    stop("\n negative weights not allowed")
  else {
    wzero <- (w == 0)
    Y.org <- Y
    X.org <- X
    offset.org <- offset
    Y <- Y * w
    X <- diag(c(w)) %*% X
    offset <- w * offset
    if (any(wzero)) {
      wpos <- !wzero
      fitted <- resid <- q1 <- q2 <- Y.org
      Y <- Y[wpos]
      X <- as.matrix(X[wpos, ])
      offset <- offset[wpos]
    }
  }
  method <- "elliptical.fit"
  elliptical.fitter <- get(method)
  offset4fit <- offset



  fit <- elliptical.fitter(X = X, Y = Y, offset = offset4fit, family = family, dispersion = dispersion, 
                             maxit = control$maxit, epsilon = control$epsilon, trace = control$trace, ...)

  if (any(wzero)) {
    nas <- is.na(fit$coef)
    fitted[wpos] <- fit$fitted.values/w[wpos]
    fitted[wzero] <- X.org[wzero, !nas] %*% as.vector(fit$coef[!nas]) + 
      if (length(offset.org) > 1) 
        offset.org[wzero]
    else 0
    fit$fitted.values <- fitted
    resid[wpos] <- fit$resid
    resid[wzero] <- (Y.org[wzero] - fitted[wzero])/sqrt(fit$dispersion)
    fit$residuals <- resid
    q1[wpos] <- fit$q1
    q2[wpos] <- fit$q2
    q1[wzero] <- family$g1(resid[wzero], df = family$df, 
                           alpha = family$alpha, mp = family$mp, epsi = family$epsi, 
                           sigmap = family$sigmap, k = family$k)
    q2[wzero] <- -2 * q1[wzero]
    fit$q1 <- q1
    fit$q2 <- q2
  }
  else fit$fitted.values <- fit$fitted.values/w
  fit$weights <- w
  names(fit$fitted.values) <- names(fit$residuals) <- names(fit$q1) <- names(fit$q2) <- NULL
  p <- dim(X)[2]
  rank <- fit$rank
  df.residuals <- length(if (exists("X.org", frame = sys.nframe())) Y.org else Y) - 
    rank - sum(w == 0) 
  asgn <- attr(if (exists("X.org", frame = sys.nframe())) X.org else X, 
               "assign")
  if (rank < p) {
    nas <- is.na(fit$coef)
    pasgn <- asgn[!nas]
    if (df.residuals > 0) 
      fit$assign.residual <- (rank + 1):length(Y)
    fit$R.assign <- pasgn
    fit$x.assign <- asgn
  }
  fit <- c(fit, list(assign = asgn, df.residuals = df.residuals, data = data,
                     family = family, user.def = user.def, formula = as.vector(attr(Terms, 
                     "formula")), terms = Terms, contrasts = attr(X, "contrasts"), 
                     control = control, call = call))
  if (y) 
    fit$y <- if (exists("Y.org", frame = sys.nframe())) 
      Y.org
  else Y
  names(fit$y) <- NULL
  if (x) 
    fit$X <- if (exists("X.org", frame = sys.nframe())) 
      X.org
  else X
  if (model) 
    fit$model <- m
  attr(fit, "class") <- c("elliptical", "glm", "lm")

  fit
}


elliptical.fit <- function (X, Y, offset, family, dispersion, 
                            maxit, epsilon, trace) 
{
  n <- nrow(X)
  if (is.null(offset)) {
    offset <- rep(0, n)
  }
  
  p <- ncol(X)
  aux.model <- glm.fit(x = X, y = Y, offset = offset, 
                       family = gaussian())
  attr(aux.model, "class") <- c("glm", "lm")
  start <- aux.model$coef
  
  
  is.null.disp <- is.null(dispersion)
  elliptical.disp <- !is.null.disp && !is.number(dispersion)
  if (is.null.disp) 
    dispersion <- (summary(aux.model)$dispersion)
  if (elliptical.disp) 
    dispersion <- (summary(aux.model)$dispersion)
  args <- resid(aux.model)/sqrt(dispersion)
  
  if (any(nas <- is.na(start))) {
    names(nas) <- dimnames(X)[[2]]
    X <- X[, !nas]
    aux.model <- glm.fit(x = X, y = Y, offset = offset, 
                         family = gaussian())
    attr(aux.model, "class") <- c("glm", "lm")
    start <- aux.model$coef
    dispersion <- (summary(aux.model)$dispersion)
  }
  
  
  iter <- 1
  error2 <- error3 <- 0
  repeat {
    if (trace) 
      cat("\n iteration", iter, ":")
    {
      w.1 <- family$g1(args, df = family$df, r = family$r, 
                       s = family$s, alpha = family$alpha, mp = family$mp, 
                       epsi = family$epsi, sigmap = family$sigmap, 
                       k = family$k)
      dg <- family$g2(args, df = family$df, r = family$r, 
                      s = family$s, alpha = family$alpha, mp = family$mp, 
                      epsi = family$epsi, sigmap = family$sigmap, 
                      k = family$k)
      fg <- family$g3(args, df = family$df, r = family$r, 
                      s = family$s, alpha = family$alpha, mp = family$mp, 
                      epsi = family$epsi, sigmap = family$sigmap, 
                      k = family$k)
      
      y.aux <- Y - offset
      w.h <- as.vector(-2 * w.1)
      aux.model <- glm.fit(x = X, y = y.aux, weights = w.h, 
                           family = gaussian())
      attr(aux.model, "class") <- c("glm", "lm")
      new.start <- coef(aux.model)
      }
    error1 <- max(abs((new.start - start)/start))
    start <- new.start
    abs.res <- Y - X %*% start - offset
    
    if (is.null.disp) {
      aux.dispersion <- dispersion
      new.dispersion <- mean((-2 * w.1) * abs.res^2)
      error2 <- abs((new.dispersion - dispersion)/dispersion)
      dispersion <- new.dispersion
    }
    old.args <- args
    args <- abs.res/sqrt(dispersion)
    if (trace) {
      loglik <- -0.5 * length(abs.res) * log((dispersion)) + 
        sum(family$g0(abs.res/sqrt(dispersion), df = family$df, 
                      s = family$s, r = family$r, alpha = family$alpha, 
                      mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                      k = family$k))
      cat(" log-likelihood =", signif(loglik, 6))
    }
    error3 <- sqrt(sum((args - old.args)^2)/max(1e-20, sum(old.args^2)))
    if ((iter == maxit) || (max(error1, error2, error3, 
                                na.rm = TRUE) < epsilon)) 
      break
    iter <- iter + 1
  }
  if (trace) 
    cat("\n")
  if (maxit > 1 && iter == maxit) 
    warning(paste("\n linear convergence not obtained in", 
                  maxit, "iterations"))
  coefs <- rep(NA, length(nas))
  coefs[!nas] <- start
  names(coefs) <- names(nas)
  names(dispersion) <- "dispersion"
  fitted <- as.vector(X %*% start + offset)
  
  
  residuals <- (Y - fitted)/sqrt(dispersion)
  w.1 <- family$g1(residuals, df = family$df, s = family$s, 
                   r = family$r, alpha = family$alpha, mp = family$mp, 
                   epsi = family$epsi, sigmap = family$sigmap, k = family$k)
  w.2 <- -2 * w.1
  if (any(w.2 < 0)) 
    cat("\n --- negative iterative weights returned! --- \n")
  if (is.null.disp) {
    rank <- dim(X)[2]
    Rnames <- dimnames(X)[[2]]
    Xd <- cbind(X, residuals)
  }
  dimnames(Xd)[[2]] <- c(Rnames, "scale")
  nn <- is.null(Rnames)
  Rnames <- list(dimnames(Xd)[[2]], dimnames(Xd)[[2]])
  R <- t(Xd) %*% Xd
  if (is.null.disp) 
    R[rank + 1, rank + 1] <- R[rank + 1, rank + 1] + length(residuals)
  attributes(R) <- list(dim = dim(R))
  if (!nn) 
    attr(R, "dimnames") <- Rnames
  loglik <- -0.5 * length(residuals) * log((dispersion)) + 
    sum(family$g0(residuals, df = family$df, s = family$s, 
                  r = family$r, alpha = family$alpha, mp = family$mp, 
                  epsi = family$epsi, sigmap = family$sigmap, k = family$k))
  names(loglik) <- NULL
  fit <- list(coefficients = coefs, dispersion = dispersion, 
              fixed = !is.null.disp, residuals = residuals, fitted.values = fitted, 
              loglik = loglik, Wg = family$g1(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), Wgder = family$g5(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), v = -2 * family$g1(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), rank = rank, R = as.matrix(R), iter = iter - 
              1, scale = 4 * family$g2(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), scaledispersion = -1 + 4 * family$g3(args, 
              df = family$df, r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), scalevariance = family$g4(args, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), df = if (charmatch(family$family, 
              "Student", F)) family$df, s = if (charmatch(family$family, 
              "Gstudent", F)) family$s, r = if (charmatch(family$family, 
              "Gstudent", F)) family$r, alpha = if (charmatch(family$family, 
              "Glogis", F)) family$alpha, mp = if (charmatch(family$family, 
              "Glogis", F)) family$m, epsi = if (charmatch(family$family, 
              "Cnormal", F)) family$epsi, sigmap = if (charmatch(family$family, 
              "Cnormal", F)) family$sigmap, k = if (charmatch(family$family, 
              "Powerexp", F)) family$k, Xmodel = matrix(Xd[, (1:rank)], 
              nrow(Xd), rank))
  fit
}

#' @rdname elliptical
#' @method print elliptical
#' @noRd
#' @export
print.elliptical <- function (x, digits = 6,...) 
{
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  coef <- x$coef
  if (any(nas <- is.na(coef))) {
    if (is.null(names(coef))) 
      names(coef) <- paste("b", 1:length(coef), sep = "")
    coef <- coef[!nas]
    cat("\nCoefficients: (", sum(nas), " not defined because of singularities)\n", 
        sep = "")
  }
  else cat("\nCoefficients:\n")
  print(coef, digits = digits, ...)
  cat("\nScale parameter: ", format(x$dispersion, digits = digits), 
      if (x$fixed) 
        " (fixed)\n"
      else "\n")
  cat("\nError distribution: ", x$family[[1]], "\n")
  rank <- x$rank
  if (is.null(rank)) 
    rank <- sum(!nas)
  nobs <- length(x$residuals) - sum(x$weights == 0)
  rdf <- x$df.resid
  if (is.null(rdf)) 
    rdf <- nobs - rank
  cat("\nDegrees of Freedom:", nobs, "Total;", rdf, "Residual\n")
  cat("-2*Log-Likelihood", format(-2 * x$loglik, digits = digits), 
      "\n")
  cat("AIC:", (2*rank - 2*x$loglik), " AICc:", (2*rank - 2*x$loglik + 2*(rank)*(rank+1)/(nobs-rank-1)), " BIC:", (log(nobs)*rank - 2*x$loglik), "\n")
  invisible(x)
}
