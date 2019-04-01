#' @title Analysis of Deviance for Elliptical Model Fits
#' @method anova elliptical
#' @description Compute an analysis of deviance table for one or more elliptical model fits.
#' @param object fit object for elliptical regression model.
#' @param dispersion the dispersion parameter for the fitting family, by default obtained from object.
#' @param test a character string representing that hypothesis test should be considered.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return An object of class \code{anova} inheriting from class \code{data.frame}
#' @references Cysneiros, F. J. A., Paula, G. A., and Galea, M. (2007). Heteroscedastic 
#' symmetrical linear models. Statistics & probability letters, 77(11), 1084-1090. 
#' \url{https://doi.org/10.1016/j.spl.2007.01.012} 
#' @seealso \code{\link{elliptical}}
#' @keywords elliptical
#' @keywords ANOVA
#' @examples
#' data(luzdat)
#' y <- luzdat$y
#' x1 <- luzdat$x1 ; x1 <- factor(x1) ; x1 <- C(x1,treatment)
#' x2 <- luzdat$x2
#' x3 <- (luzdat$x2)^2
#' luz <- data.frame(y,x1,x2,x3)
#' elliptical.fitt <- elliptical(y ~ x1+x2+x3, family = Student(df=5),
#' data=luz)
#' anova(elliptical.fitt)
#' @export
 
anova.elliptical <- function(object, dispersion = NULL, test = c("Chisq"),...) 
{
  test <- match.arg(test)
  margs <- function(...) nargs()
  if (margs(...)) 
    return(anova.ellipticallist(list(object, ...), test = test))
  Terms <- object$terms
  term.labels <- attr(Terms, "term.labels")
  nt <- length(term.labels)
  m <- model.frame(object)
  family <- family(object)[[1]]
  y <- model.extract(m, "response")
  loglik <- double(nt + 1)
  df.res <- loglik
  if (nt) {
    loglik[nt + 1] <- -2 * object$loglik
    df.res[nt + 1] <- object$df.residual
    fit <- object
    for (iterm in seq(from = nt, to = 1, by = -1)) {
      ww <- fit$weights
      argslist <- list(object = fit, formula = eval(parse(text = paste("~ . -", 
                                                                       term.labels[iterm]))))
      fit <- do.call("update", argslist)
      loglik[iterm] <- -2 * fit$loglik
      df.res[iterm] <- fit$df.residual
    }
    dev <- c(NA, -diff(loglik))
    df <- c(NA, -diff(df.res))
  }
  else {
    loglik[1] <- -2 * object$loglik
    df.res[1] <- dim(y)[1] - attr(Terms, "intercept")
    dev <- df <- as.numeric(NA)
  }
  heading <- c("Analysis of Deviance Table\n", paste("Error distribution:", 
                                                     family), paste("Response:", as.character(formula(object))[2], 
                                                                    "\n", sep = ""), "Terms added sequentially (first to last)")
  if (is.null(dispersion)) {
    dispersion <- 1
    df.dispersion <- if (dispersion == 1) 
      Inf
    else object$df.residual
  }
  else df.scale <- Inf
  aod <- data.frame(Df = df, Deviance = dev, `Resid. Df` = df.res, 
                    `-2*LL` = loglik, row.names = c("NULL", term.labels), 
                    check.names = F)
  attr(aod, "heading") <- heading
  class(aod) <- c("anova", "data.frame")
  if (is.null(test)) 
    return(aod)
  else aod <- stat.anova(aod, test, scale = 1, df.scale = df.dispersion, 
                         n = nrow(y))
  structure(aod, heading = heading, class = c("anova", "data.frame"))
}


anova.ellipticallist <- function (object, test = c("Chisq")) 
{
  diff.term <- function(term.labels, i) {
    t1 <- term.labels[[1]]
    t2 <- term.labels[[2]]
    m1 <- match(t1, t2, F)
    m2 <- match(t2, t1, F)
    if (all(m1)) {
      if (all(m2)) 
        return("=")
      else return(paste(c("", t2[-m1]), collapse = "+"))
    }
    else {
      if (all(m2)) 
        return(paste(c("", t1[-m2]), collapse = "-"))
      else return(paste(i - 1, i, sep = " vs. "))
    }
  }
  test <- match.arg(test)
  rt <- length(object)
  if (rt == 1) {
    object <- object[[1]]
    UseMethod("anova")
  }
  forms <- sapply(object, function(x) as.character(formula(x)))
  subs <- as.logical(match(forms[2, ], forms[2, 1], F))
  if (!all(subs)) 
    warning("Some fit objects deleted because response differs from the first model")
  if (sum(subs) == 1) 
    stop("The first model has a different response from the rest")
  forms <- forms[, subs]
  object <- object[subs]
  dfres <- sapply(object, "[[", "df.resid")
  m2loglik <- -2 * sapply(object, "[[", "loglik")
  tl <- lapply(object, labels)
  rt <- length(m2loglik)
  effects <- character(rt)
  for (i in 2:rt) effects[i] <- diff.term(tl[c(i - 1, i)], 
                                          i)
  dm2loglik <- -diff(m2loglik)
  ddf <- -diff(dfres)
  family <- family(object[[1]])[[1]]
  fixed <- object[[1]]$fixed
  heading <- c("Analysis of Deviance Table\n", paste("Error distribution:", 
                                                     family), paste("Response:", forms[2, 1], "\n"))
  nmodels <- length(object)
  #  varirest <- lapply(object, function(x) paste(deparse(x$call$restrict), 
  #                                               collapse = "\n"))
  variables <- lapply(object, function(x) paste(deparse(formula(x)), 
                                                collapse = "\n"))
  #  topnote <- paste("Model ", format(1:nmodels), ": ", variables, 
  #                   "     restrict =", varirest, sep = "", collapse = "\n")
  aod <- data.frame(Terms = forms[3, ], `Resid. Df` = dfres, 
                    `-2*LL` = m2loglik, Test = effects, Df = c(NA, abs(ddf)), 
                    Deviance = c(NA, abs(dm2loglik)), check.names = F)
  attr(aod, "heading") <- heading
  if (!is.null(test)) {
    n <- length(object[[1]]$residuals)
    o <- order(dfres)
    aod <- stat.anova(aod, test, 1, dfres[o[1]], n)
  }
  else aod
  #  structure(aod, heading = c(heading, topnote), class = c("anova",
  structure(aod, heading = heading, class = c("anova",
                                              "data.frame"))
}


