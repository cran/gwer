#' @title Diagnostic Plots for Elliptical Regression Models
#' @import utils
#' @import graphics
#' @description This function produces diagnostic measures plots for elliptical regression models.
#' @param ellipticalfit fit object for elliptical regression model.
#' @param ellipticaldiag objects containing the diagnostic measures, by default obtained from object.
#' @param weighting type of model weighting used.
#' @param which an optional numerical that indicates which plot is returned.
#' @param subset optional vector specifying a subset of observations to be used in the fitting process.
#' @param iden a logical value used to identify observations. if TRUE the observations can be identified in the graphic window.
#' @param labels a optinal vector specifying a labels plots.
#' @param ret a logical value to indicate funtion returns. If TRUE the return of the function will be to the diagnostic measures used.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return If \code{ret} is true, returns a list of diagnostic arrays (see \code{elliptical.diag} for more details).
#' @references Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
#' symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
#' \url{https://doi.org/10.1016/j.spl.2005.04.033}
#' @seealso \code{\link{elliptical}}, \code{\link{elliptical.diag}}
#' @keywords elliptical
#' @keywords diagnostic methods
#' @examples
#' data(luzdat)
#' y <- luzdat$y
#' x1 <- luzdat$x1 ; x1 <- factor(x1) ; x1 <- C(x1,treatment)
#' x2 <- luzdat$x2
#' x3 <- (luzdat$x2)^2
#' luz <- data.frame(y,x1,x2,x3)
#' elliptical.fitt <- elliptical(y ~ x1+x2+x3, family = Student(df=5),
#' data=luz)
#' elliptical.diag.plots(elliptical.fitt, which=3)
#' @export


elliptical.diag.plots <- function (ellipticalfit, ellipticaldiag = NULL, weighting, which, 
          subset = NULL, iden = F, labels = NULL, ret = F,...) 
{
  if (is.null(ellipticaldiag)) {
    if (missing(weighting)) {
      family <- ellipticalfit$family
      user.def <- ellipticalfit$user.def
      f.name <- family[[1]]
      weighting <- "observed"
    }
    ellipticaldiag <- elliptical.diag(ellipticalfit, weighting = weighting)
  }
  if (is.null(subset)) 
    subset <- c(1:length(ellipticaldiag$GL))
  else if (is.logical(subset)) 
    subset <- (1:length(subset))[subset]
  else if (is.numeric(subset) && all(subset < 0)) 
    subset <- (1:(length(subset) + length(ellipticaldiag$GL)))[subset]
  else if (is.character(subset)) {
    if (is.null(labels)) 
      labels <- subset
    subset <- seq(along = subset)
  }
  w <- if (is.null(ellipticalfit$weights)) 
    rep(1, length(ellipticalfit$residuals))
  else ellipticalfit$weights
  wzero <- (w == 0)
  choices <- c("All", "Response residual against fitted values", 
               "Response residual against index", "Standardized residual against fitted values", 
               "Standardized residual against index", "QQ-plot of  response residuals", 
               "QQ-plot of  Standardized residuals", "Generalized Leverage ", 
               "Ci against index ", "|Lmax| against index (local influence on coefficients)", 
               "Bii against index\n")
  tmenu <- paste("plot:", choices)
  if (missing(which)) 
    pick <- menu(tmenu, title = "\n Make a plot selection (or 0 to exit)\n")
  else if (!match(which, 2:11, nomatch = F)) 
    stop("choice not valid")
  else pick <- which
  if (pick == 0) 
    stop(" no graph required ! ")
  repeat {
    switch(pick, `1` = {
      par(pty = "s", mfrow = c(1, 1))
      close.screen(all.screens = T)
      split.screen(c(2, 2))
      screen(1)
      x1 <- ellipticalfit$fitted.values[!wzero]
      y1 <- ellipticaldiag$resid
      plot(x1, y1, xlab = "Fitted values", ylab = "Response residual", 
           ...)
      screen(2)
      x2 <- 1:length(ellipticalfit$fitted.values[!wzero])
      y2 <- ellipticaldiag$resid
      plot(x2, y2, xlab = "Index", ylab = "Response residual", 
           ...)
      screen(3)
      x3 <- ellipticalfit$fitted.values[!wzero]
      y3 <- ellipticaldiag$rs
      plot(x3, y3, xlab = "Fitted values", ylab = "Standardized residual", 
           ...)
      screen(4)
      x4 <- 1:length(ellipticalfit$fitted.values[!wzero])
      y4 <- ellipticaldiag$rs
      plot(x4, y4, xlab = "Index", ylab = "Standardized residual", 
           ...)
      xx <- list(x1, x2, x3, x4)
      yy <- list(y1, y2, y3, y4)
      if (is.null(labels)) labels <- names(model.extract(model.frame(ellipticalfit), 
                                                         "response"))
      yes <- iden
      while (yes) {
        cat("****************************************************\n")
        cat("Please Input a screen number  (1,2,3, or 4)\n")
        cat("0 will terminate the function \n")
        num <- scan(n = 1)
        if ((length(num) > 0) && ((num == 1) || (num == 
                                                 2) || (num == 3) || (num == 4))) {
          cat(paste("Interactive Identification for screen", 
                    num, "\n"))
          cat("left button = Identify, center button = Exit\n")
          screen(num, new = F)
          identify(xx[[num]], yy[[num]], labels, ...)
        } else yes <- F
      }
      close.screen(all.screens = T)
      par(ask = T)
      split.screen(figs = c(2, 2))
      screen(1)
      y5 <- ellipticaldiag$resid
      x5 <- qnorm(ppoints(length(y5)))
      x5 <- x5[rank(y5)]
      .lim <- c(min(x5, y5), max(x5, y5))
      plot(x5, y5, xlab = paste("Quantiles of standard normal"), 
           ylab = "Ordered response residual", xlim = .lim, 
           ylim = .lim, ...)
      abline(0, 1, lty = 2)
      screen(2)
      y6 <- ellipticaldiag$rs
      x6 <- qnorm(ppoints(length(y6)))
      x6 <- x6[rank(y6)]
      .lim <- c(min(x6, y6), max(x6, y6))
      plot(x6, y6, xlab = paste("Quantiles of standard normal"), 
           ylab = "Ordered standardized residual", xlim = .lim, 
           ylim = .lim, ...)
      abline(0, 1, lty = 2)
      screen(3)
      y7 <- ellipticaldiag$GL
      x7 <- 1:length(ellipticalfit$fitted.values[!wzero])
      plot(x7, y7, xlab = "Index", ylab = "Generalized leverage ", 
           ...)
      screen(4)
      y8 <- ellipticaldiag$Ci
      x8 <- 1:length(ellipticalfit$fitted.values[!wzero])
      plot(x8, y8, xlab = "Index", ylab = "Ci", ...)
      xx <- list(x5, x6, x7, x8)
      yy <- list(y5, y6, y7, y8)
      if (is.null(labels)) labels <- names(model.extract(model.frame(ellipticalfit), 
                                                         "response"))
      yes <- iden
      while (yes) {
        cat("****************************************************\n")
        cat("Please Input a screen number  (1,2,3, or 4)\n")
        cat("0 will terminate the function \n")
        num <- scan(n = 1)
        if ((length(num) > 0) && ((num == 1) || (num == 
                                                 2) || (num == 3) || (num == 4))) {
          cat(paste("Interactive Identification for screen", 
                    num, "\n"))
          cat("left button = Identify, center button = Exit\n")
          screen(num, new = F)
          identify(xx[[num]], yy[[num]], labels, ...)
        } else yes <- F
      }
      close.screen(all.screens = T)
      par(ask = T)
      split.screen(figs = c(2, 2))
      screen(1)
      y9 <- ellipticaldiag$dmax
      x9 <- 1:length(ellipticalfit$fitted.values[!wzero])
      plot(x9, y9, xlab = "Index", ylab = "|Lmax|", ...)
      screen(2)
      y10 <- ellipticaldiag$Bi
      x10 <- 1:length(ellipticalfit$fitted.values[!wzero])
      plot(x10, y10, xlab = "Index", ylab = "Bii", ...)
      xx <- list(x9, x10)
      yy <- list(y9, y10)
      if (is.null(labels)) labels <- names(model.extract(model.frame(ellipticalfit), 
                                                         "response"))
      yes <- iden
      while (yes) {
        cat("****************************************************\n")
        cat("Please Input a screen number (1, 2 , 3 or 4)\n")
        cat("0 will terminate the function \n")
        num <- scan(n = 1)
        if ((length(num) > 0) && ((num == 1) || (num == 
                                                 2) || (num == 3) || (num == 4))) {
          cat(paste("Interactive Identification for screen", 
                    num, "\n"))
          cat("left button = Identify, center button = Exit\n")
          screen(num, new = F)
          identify(xx[[num]], yy[[num]], labels, ...)
        } else yes <- F
      }
      close.screen(all.screens = T)
      par(ask = F)
    }, `2` = {
      par(pty = "s")
      x2 <- ellipticalfit$fitted.values[!wzero]
      y2 <- ellipticaldiag$resid
      plot(x2, y2, xlab = "Fitted values", ylab = "Response residual", 
           ...)
      xx <- list(x2)
      yy <- list(y2)
    }, `3` = {
      par(pty = "s")
      x3 <- 1:length(ellipticalfit$fitted.values[!wzero])
      y3 <- ellipticaldiag$resid
      plot(x3, y3, xlab = "Index", ylab = "Response residual", 
           ...)
      xx <- list(x3)
      yy <- list(y3)
    }, `4` = {
      par(pty = "s")
      x4 <- ellipticalfit$fitted.values[!wzero]
      y4 <- ellipticaldiag$rs
      plot(x4, y4, xlab = "Fitted values", ylab = "Standardized residual", 
           ...)
      xx <- list(x4)
      yy <- list(y4)
    }, `5` = {
      par(pty = "s")
      x5 <- 1:length(ellipticalfit$fitted.values[!wzero])
      y5 <- ellipticaldiag$rs
      plot(x5, y5, xlab = "Index", ylab = "Standardized residual", 
           ...)
      xx <- list(x5)
      yy <- list(y5)
    }, `6` = {
      par(pty = "s")
      y6 <- ellipticaldiag$resid
      x6 <- qnorm(ppoints(length(y6)))
      x6 <- x6[rank(y6)]
      .lim <- c(min(x6, y6), max(x6, y6))
      plot(x6, y6, xlab = "Quantiles of standard normal", 
           ylab = "Ordered response residual", xlim = .lim, 
           ylim = .lim, ...)
      abline(0, 1, lty = 2)
      xx <- list(x6)
      yy <- list(y6)
    }, `7` = {
      par(pty = "s")
      y7 <- ellipticaldiag$rs
      x7 <- qnorm(ppoints(length(y7)))[rank(y7)]
      .lim <- c(min(x7, y7), max(x7, y7))
      plot(x7, y7, xlab = paste("Quantiles of standard normal"), 
           ylab = "Ordered standardized residual", xlim = .lim, 
           ylim = .lim, ...)
      abline(0, 1, lty = 2)
      xx <- list(x7)
      yy <- list(y7)
    }, `8` = {
      par(pty = "s")
      y8 <- ellipticaldiag$GL
      x8 <- 1:length(ellipticalfit$fitted.values[!wzero])
      plot(x8, y8, xlab = "Index", ylab = "Generalized leverage ", 
           ...)
      xx <- list(x8)
      yy <- list(y8)
    }, `9` = {
      par(pty = "s")
      y9 <- ellipticaldiag$Ci
      x9 <- 1:length(ellipticalfit$fitted.values[!wzero])
      plot(x9, y9, xlab = "Index", ylab = "Ci", ...)
      xx <- list(x9)
      yy <- list(y9)
    }, `10` = {
      par(pty = "s")
      y10 <- ellipticaldiag$dmax
      x10 <- 1:length(ellipticalfit$fitted.values[!wzero])
      plot(x10, y10, xlab = "Index", ylab = "|Lmax|", 
           ...)
      xx <- list(x10)
      yy <- list(y10)
    }, `11` = {
      par(pty = "s")
      y11 <- ellipticaldiag$Bi
      x11 <- 1:length(ellipticalfit$fitted.values[!wzero])
      plot(x11, y11, xlab = "Index", ylab = "Bii", ...)
      xx <- list(x11)
      yy <- list(y11)
    })
    if ((!(pick == 1))) {
      if (is.null(labels)) 
        labels <- names(model.extract(model.frame(ellipticalfit), 
                                      "response"))
      yes <- iden
      while (yes) {
        cat("****************************************************\n")
        cat("Interactive Identification\n")
        cat("left button = Identify, center button = Exit\n")
        identify(xx[[1]], yy[[1]], labels, ...)
        yes <- F
      }
    }
    if (missing(which)) 
      pick <- menu(tmenu, title = "\n Make a plot selection (or 0 to exit)\n")
    if ((pick == 0) || !missing(which)) {
      invisible(close.screen(all.screens = T))
      break
    }
  }
  if (ret) 
    ellipticaldiag
  else invisible()
}
