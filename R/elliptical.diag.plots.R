#' @title Diagnostic Plots for Elliptical Regression Models
#' @import utils
#' @import graphics
#' @description This function generate diagnostic measures plots for the fitted elliptical regression models.
#' @param object an object with the result of the fitted elliptical regression model.
#' @param ellipticaldiag object list containing the diagnostic measures. By default it is obtained from the object, but can be calculated using \code{\link{elliptical.diag}}.
#' @param which an optional numeric value with the number of only plot that must be returned.
#' @param subset an optional numeric vector specifying a subset of observations to be used in the fitting process.
#' @param iden a logical value used to identify observations. If \code{TRUE} the observations are identified by user in the graphic window.
#' @param labels a optional string vector specifying a labels plots.
#' @param ret a logical value used to return the diagnostic measures computing. If \code{TRUE} the diagnostic measures are returned (see \code{\link{elliptical.diag}} for more details).
#' @param ... graphics parameters to be passed to the plotting routines.
#' @return Return an interactive menu with eleven options to make plots. This menu contains the follows graphics:
#' 1: plot: All.
#' 2: plot: Response residual against fitted values.
#' 3: plot: Response residual against index.
#' 4: plot: Standardized residual against fitted values.
#' 5: plot: Standardized residual against index.
#' 6: plot: QQ-plot of response residuals.
#' 7: plot: QQ-plot of standardized residuals.
#' 8: plot: Generalized leverage.
#' 9: plot: Total local influence index plot for response perturbation.
#' 10: plot: Total local influence index plot scale perturbation.
#' 11: plot: Total local influence index plot case-weight perturbation.
#' @references Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
#' symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
#' \doi{10.1016/j.spl.2005.04.033}
#' @seealso \code{\link{elliptical}}, \code{\link{elliptical.diag}}
#' @keywords Elliptical regression models
#' @keywords Diagnostic methods
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

elliptical.diag.plots <- function (object, ellipticaldiag = NULL, which, 
          subset = NULL, iden = FALSE, labels = NULL, ret = FALSE, ...) 
{
  if (is.null(ellipticaldiag)) {
    family <- object$family
    user.def <- object$user.def
    f.name <- family[[1]]
    ellipticaldiag <- elliptical.diag(object)
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
  nobs = length(object$residuals)
  w <- if (is.null(object$weights)) 
    rep(1, nobs)
  else object$weights
  wzero <- (w == 0)
  choices <- c("All", "Response residual against fitted values", 
               "Response residual against index", "Standardized residual against fitted values", 
               "Standardized residual against index", "QQ-plot of response residuals", 
               "QQ-plot of standardized residuals", "Generalized leverage", 
               "Total local influence index plot for response perturbation", 
               "Total local influence index plot scale perturbation", 
               "Total local influence index plot case-weight perturbation\n")
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
      x1 <- object$fitted.values[!wzero]
      y1 <- ellipticaldiag$rr
      plot(x1, y1, xlab = "Fitted values", ylab = "Response residual", 
           ...)
      screen(2)
      x2 <- 1:length(object$fitted.values[!wzero])
      y2 <- ellipticaldiag$rr
      plot(x2, y2, xlab = "Index", ylab = "Response residual", 
           ...)
      screen(3)
      x3 <- object$fitted.values[!wzero]
      y3 <- ellipticaldiag$rs
      plot(x3, y3, xlab = "Fitted values", ylab = "Standardized residual", 
           ...)
      screen(4)
      x4 <- 1:length(object$fitted.values[!wzero])
      y4 <- ellipticaldiag$rs
      plot(x4, y4, xlab = "Index", ylab = "Standardized residual", 
           ...)
      xx <- list(x1, x2, x3, x4)
      yy <- list(y1, y2, y3, y4)
      if (is.null(labels)) labels <- names(model.extract(model.frame(object), 
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
      y5 <- ellipticaldiag$rr
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
      x7 <- 1:length(object$fitted.values[!wzero])
      plot(x7, y7, xlab = "Index", ylab = "Generalized leverage ", ylim = c(0,1), 
           ...) ; abline(2/nobs,0)
      screen(4)
      y8 <- ellipticaldiag$Lmaxr
      x8 <- 1:length(object$fitted.values[!wzero])
      plot(x8, y8, xlab = "Index", ylab = "|Lmax|", ylim = c(0,1), ...) ; abline(2/nobs,0)
      xx <- list(x5, x6, x7, x8)
      yy <- list(y5, y6, y7, y8)
      if (is.null(labels)) labels <- names(model.extract(model.frame(object), 
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
      y9 <- ellipticaldiag$Cic
      x9 <- 1:length(object$fitted.values[!wzero])
      plot(x9, y9, xlab = "Index", ylab = "Ci", ylim = c(0,1), ...) ; abline(2/nobs,0)
      screen(2)
      y10 <- ellipticaldiag$Cih
      x10 <- 1:length(object$fitted.values[!wzero])
      plot(x10, y10, xlab = "Index", ylab = "Ci", ylim = c(0,1), ...) ; abline(2/nobs,0)
      xx <- list(x9, x10)
      yy <- list(y9, y10)
      if (is.null(labels)) labels <- names(model.extract(model.frame(object), 
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
      x2 <- object$fitted.values[!wzero]
      y2 <- ellipticaldiag$rr
      plot(x2, y2, xlab = "Fitted values", ylab = "Response residual", 
           ...)
      xx <- list(x2)
      yy <- list(y2)
    }, `3` = {
      par(pty = "s")
      x3 <- 1:length(object$fitted.values[!wzero])
      y3 <- ellipticaldiag$rr
      plot(x3, y3, xlab = "Index", ylab = "Response residual", 
           ...)
      xx <- list(x3)
      yy <- list(y3)
    }, `4` = {
      par(pty = "s")
      x4 <- object$fitted.values[!wzero]
      y4 <- ellipticaldiag$rs
      plot(x4, y4, xlab = "Fitted values", ylab = "Standardized residual", 
           ...)
      xx <- list(x4)
      yy <- list(y4)
    }, `5` = {
      par(pty = "s")
      x5 <- 1:length(object$fitted.values[!wzero])
      y5 <- ellipticaldiag$rs
      plot(x5, y5, xlab = "Index", ylab = "Standardized residual", 
           ...)
      xx <- list(x5)
      yy <- list(y5)
    }, `6` = {
      par(pty = "s")
      y6 <- ellipticaldiag$rr
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
      x8 <- 1:length(object$fitted.values[!wzero])
      plot(x8, y8, xlab = "Index", ylab = "Generalized leverage ", 
           ylim = c(0,1), ...) ; abline(2/nobs,0)
      xx <- list(x8)
      yy <- list(y8)
    }, `9` = {
      par(pty = "s")
      y9 <- ellipticaldiag$Lmaxr
      x9 <- 1:length(object$fitted.values[!wzero])
      plot(x9, y9, xlab = "Index", ylab = "|Lmax|", ylim = c(0,1),
           ...) ; abline(2/nobs,0)
      xx <- list(x9)
      yy <- list(y9)
    }, `10` = {
      par(pty = "s")
      y10 <- ellipticaldiag$Cic
      x10 <- 1:length(object$fitted.values[!wzero])
      plot(x10, y10, xlab = "Index", ylab = "Ci", 
           ylim = c(0,1), ...) ; abline(2/nobs,0)
      xx <- list(x10)
      yy <- list(y10)
    }, `11` = {
      par(pty = "s")
      y11 <- ellipticaldiag$Cih
      x11 <- 1:length(object$fitted.values[!wzero])
      plot(x11, y11, xlab = "Index", ylab = "Ci", ylim = c(0,1),
           ...) ; abline(2/nobs,0)
      xx <- list(x11)
      yy <- list(y11)
    })
    if ((!(pick == 1))) {
      if (is.null(labels)) 
        labels <- names(model.extract(model.frame(object), 
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
