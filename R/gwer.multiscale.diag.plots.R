#' @title Diagnostic Plots for Multiscale Geographically Weighted Elliptical Regression Models
#' @import utils
#' @import graphics
#' @description This function generate diagnostic measures plots for the fitted multiscale geographically weighted elliptical regression models.
#' @param object an object with the result of the fitted multiscale geographically weighted elliptical regression models.
#' @param mgwerdiag object list containing the diagnostic measures. By default it is obtained from the object, but can be calculated using \code{\link{gwer.multiscale.diag}}.
#' @param which an optional numeric value with the number of only plot that must be returned.
#' @param subset an optional numeric vector specifying a subset of observations to be used in the fitting process.
#' @param iden a logical value used to identify observations. If \code{TRUE} the observations are identified by user in the graphic window.
#' @param labels a optional string vector specifying a labels plots.
#' @param ret a logical value used to return the diagnostic measures computing. If \code{TRUE} the diagnostic measures are returned (see \code{\link{gwer.multiscale.diag}} for more details).
#' @param ... graphics parameters to be passed to the plotting routines.
#' @return Return an interactive menu with eleven options to make plots. This menu contains the follows graphics:
#' 1: plot: All.
#' 2: plot: Response residual against fitted values.
#' 3: plot: Response residual against index.
#' 4: plot: Quantile residual against fitted values.
#' 5: plot: Quantile residual against index.
#' 6: plot: QQ-plot of response residuals.
#' 7: plot: QQ-plot of Quantile residuals.
#' 8: plot: Generalized leverage.
#' 9: plot: Total local influence index plot for response perturbation.
#' 10: plot: Total local influence index plot scale perturbation.
#' 11: plot: Total local influence index plot case-weight perturbation.
#' @references Brunsdon, C., Fotheringham, A. S. and Charlton, M. E. (1996). 
#' Geographically weighted regression: a method for exploring spatial nonstationarity.
#' Geographical analysis, 28(4), 281-298. \doi{10.1111/j.1538-4632.1996.tb00936.x}
#' @references Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
#' symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
#' \doi{10.1016/j.spl.2005.04.033}
#' @seealso \code{\link{gwer.multiscale}}, \code{\link{gwer.multiscale.diag}}
#' @keywords Multiscale Geographically Weighted Regression
#' @keywords Elliptical models
#' @keywords Diagnostic methods
#' @examples
#' \donttest{
#' data(georgia, package = "spgwr")
#' fit.formula <- PctBach ~ TotPop90 + PctRural + PctFB + PctPov
#' gwer.bw.t <- bw.gwer(fit.formula, data = gSRDF, family = Student(3), adapt = TRUE)
#' msgwr.fit.t <- gwer.multiscale(fit.formula, family = Student(3), data = gSRDF, 
#'                                bws0 = rep(gwer.bw.t, 5), hatmatrix = TRUE, 
#'                                adaptive = TRUE)
#' gwer.multiscale.diag.plots(msgwr.fit.t, which=3)
#' }
#' @export


gwer.multiscale.diag.plots <- function (object, mgwerdiag = NULL, which, subset = NULL, 
                              iden = F, labels = NULL, ret = F,...) 
{
  if(!object$hatmatrix) 
    stop("Diagnostic measures not applicable - regression points different from observed points")
  
  if (is.null(mgwerdiag)) {
    family <- object$family
    user.def <- object$lm$user.def
    f.name <- family[[1]]
    mgwerdiag <- gwer.multiscale.diag(object)
  }
  if (is.null(subset)) 
    subset <- c(1:length(mgwerdiag$H))
  else if (is.logical(subset)) 
    subset <- (1:length(subset))[subset]
  else if (is.numeric(subset) && all(subset < 0)) 
    subset <- (1:(length(subset) + length(mgwerdiag$H)))[subset]
  else if (is.character(subset)) {
    if (is.null(labels)) 
      labels <- subset
    subset <- seq(along = subset)
  }
  nobs = length(object$lm$residuals)
  w <- if (is.null(object$lm$weights)) 
    rep(1, nobs)
  else object$lm$weights
  wzero <- (w == 0)
  p <- object$lm$rank 
  
  choices <- c("All", "Response residual against fitted values", 
               "Response residual against index", "Quantile residual against fitted values", 
               "Quantile residual against index", "QQ-plot of response residuals", 
               "QQ-plot of quantile residuals", "Generalized leverage", 
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
      split.screen(c(2, 4))
      screen(1)
      x1 <- object$fitted[!wzero]
      y1 <- mgwerdiag$rr
      plot(x1, y1, xlab = "Fitted values", ylab = "Response residual", 
           ...)
      screen(2)
      coord <- coordinates(object$SDF)
      distcoord <- knn2nb(knearneigh(coord, longlat = object$longlat))
      col.test <- nb2listw(distcoord, style="W")
      y2 <- mgwerdiag$rr
      moran.plot(y2, col.test, xlab = "Response residual", ylab = "Spatially lagged response residual", ...)
      screen(3)
      y3 <- mgwerdiag$rr
      x3 <- qnorm(ppoints(length(y3)))
      x3 <- x3[rank(y3)]
      .lim <- c(min(x3, y3), max(x3, y3))
      plot(x3, y3, xlab = paste("Quantiles of standard normal"), 
           ylab = "Ordered response residual", xlim = .lim, 
           ylim = .lim, ...)
      abline(0, 1, lty = 2)
      screen(4)
      x4 <- object$fitted[!wzero]
      y4 <- mgwerdiag$rq
      plot(x4, y4, xlab = "Fitted values", ylab = "Quatile residual", 
           ...)
      screen(5)
      coord <- coordinates(object$SDF)
      distcoord <- knn2nb(knearneigh(coord, longlat = object$longlat))
      col.test <- nb2listw(distcoord, style="W")
      y5 <- mgwerdiag$rq
      moran.plot(y5, col.test, xlab = "Quatile  residual", ylab = "Spatially lagged quatile  residual", ...)
      screen(6)
      y6 <- mgwerdiag$rq
      x6 <- qnorm(ppoints(length(y6)))
      x6 <- x6[rank(y6)]
      .lim <- c(min(x6, y6), max(x6, y6))
      plot(x6, y6, xlab = paste("Quantiles of standard normal"), 
           ylab = "Ordered Quatile residual", xlim = .lim, 
           ylim = .lim, ...)
      abline(0, 1, lty = 2)
      screen(7)
      y7 <- mgwerdiag$h/sum(mgwerdiag$h)
      x7 <- 1:length(object$fitted[!wzero])
      plot(x7, y7, xlab = "Index", ylab = "Leverage", ylim = c(0,1), 
           ...) ; abline(2/nobs,0)
      
      xx <- list(x1, x3, x4, x6, x7)
      yy <- list(y1, y3, y4, y6, y7)
      if (is.null(labels)) labels <- names(object$lm$y)
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
      split.screen(figs = c(2, ceiling(p/2)))
      x8 <- y8 <- matrix(0, nobs, p)
      for(k in 1:p){
        screen(k)
        y8[, k] <- mgwerdiag$Cic[[k]]/sum(mgwerdiag$Cic[[k]])
        x8[, k] <- 1:length(object$fitted[!wzero])
        plot(x8[, k], y8[, k], xlab = "Index", ylab = "Cic", ylim = c(0,1), ...) ; abline(2/nobs,0)
      }
      xx <- apply(x8, 2, list)#list(x5, x6, x7, x8)
      yy <- apply(y8, 2, list)#list(y5, y6, y7, y8)
      if (is.null(labels)) labels <- names(object$lm$y)
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
      split.screen(figs = c(2, ceiling(p/2)))
      x9 <- y9 <- matrix(0, nobs, p)
      for(k in 1:p){
        screen(k)
        y9[, k] <- mgwerdiag$Cih[[k]]/sum(mgwerdiag$Cih[[k]])
        x9[, k] <- 1:length(object$fitted[!wzero])
        plot(x9[, k], y9[, k], xlab = "Index", ylab = "Cih", ylim = c(0,1), ...) ; abline(2/nobs,0)
      }
      xx <- apply(x9, 2, list)#list(x5, x6, x7, x8)
      yy <- apply(y9, 2, list)#list(y5, y6, y7, y8)
      if (is.null(labels)) labels <- names(object$lm$y)
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
      par(ask = T)
      split.screen(figs = c(2, ceiling(p/2)))
      x10 <- y10 <- matrix(0, nobs, p)
      for(k in 1:p){
        screen(k)
        y10[, k] <- mgwerdiag$Lmaxr[[k]]/sum(mgwerdiag$Lmaxr[[k]])
        x10[, k] <- 1:length(object$fitted[!wzero])
        plot(x10[, k], y10[, k], xlab = "Index", ylab = "|Lmax|", ylim = c(0,1), ...) ; abline(2/nobs,0)
      }
      xx <- apply(x10, 2, list)#list(x5, x6, x7, x8)
      yy <- apply(y10, 2, list)#list(y5, y6, y7, y8)
      if (is.null(labels)) labels <- names(object$lm$y)
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
      x2 <- object$fitted[!wzero]
      y2 <- mgwerdiag$rr
      plot(x2, y2, xlab = "Fitted values", ylab = "Response residual", 
           ...)
      xx <- list(x2)
      yy <- list(y2)
    }, `3` = {
      par(pty = "s")
      coord <- coordinates(object$SDF)
      distcoord <- knn2nb(knearneigh(coord, longlat = F))
      col.test <- nb2listw(distcoord, style="W")
      y3 <- mgwerdiag$rr
      moran.plot(y3, col.test, xlab = "Response residual", ylab = "Spatially lagged response residual", ...)
    }, `4` = {
      par(pty = "s")
      y4 <- mgwerdiag$rr
      x4 <- qnorm(ppoints(length(y4)))
      x4 <- x4[rank(y4)]
      .lim <- c(min(x4, y4), max(x4, y4))
      plot(x4, y4, xlab = "Quantiles of standard normal", 
           ylab = "Ordered response residual", xlim = .lim, 
           ylim = .lim, ...)
      abline(0, 1, lty = 2)
      xx <- list(x4)
      yy <- list(y4)
    }, `5` = {
      par(pty = "s")
      x5 <- object$fitted[!wzero]
      y5 <- mgwerdiag$rq
      plot(x5, y5, xlab = "Fitted values", ylab = "Standardized residual", 
           ...)
      xx <- list(x5)
      yy <- list(y5)
    }, `6` = {
      par(pty = "s")
      coord <- coordinates(object$SDF)
      distcoord <- knn2nb(knearneigh(coord, longlat = object$longlat))
      col.test <- nb2listw(distcoord, style="W")
      y6 <- mgwerdiag$rr
      moran.plot(y6, col.test, xlab = "Standardized residual", ylab = "Spatially lagged standardized residual", ...)
    }, `7` = {
      par(pty = "s")
      y7 <- mgwerdiag$rq
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
      y8 <- mgwerdiag$h/sum(mgwerdiag$h)
      x8 <- 1:length(object$fitted[!wzero])
      plot(x8, y8, xlab = "Index", ylab = "Generalized leverage ", ylim = c(0,1), 
           ...) ; abline(2/nobs,0)
      xx <- list(x8)
      yy <- list(y8)
    }, `9` = {
      par(pty = "s")
      split.screen(figs = c(2, ceiling(p/2)))
      x9 <- y9 <- matrix(0, nobs, p)
      for(k in 1:p){
        screen(k)
        y9[, k] <- mgwerdiag$Cic[[k]]/sum(mgwerdiag$Cic[[k]])
        x9[, k] <- 1:length(object$fitted[!wzero])
        plot(x9[, k], y9[, k], xlab = "Index", ylab = "Cic", ylim = c(0,1), ...) ; abline(2/nobs,0)
      }
      xx <- apply(x9, 2, list)#list(x5, x6, x7, x8)
      yy <- apply(y9, 2, list)#list(y5, y6, y7, y8)
    }, `10` = {
      par(pty = "s")
      split.screen(figs = c(2, ceiling(p/2)))
      x10 <- y10 <- matrix(0, nobs, p)
      for(k in 1:p){
        screen(k)
        y10[, k] <- mgwerdiag$Cih[[k]]/sum(mgwerdiag$Cih[[k]])
        x10[, k] <- 1:length(object$fitted[!wzero])
        plot(x10[, k], y10[, k], xlab = "Index", ylab = "Cih", ylim = c(0,1), ...) ; abline(2/nobs,0)
      }
      xx <- apply(x10, 2, list)#list(x5, x6, x7, x8)
      yy <- apply(y10, 2, list)#list(y5, y6, y7, y8)
    }, `11` = {
      par(pty = "s")
      split.screen(figs = c(2, ceiling(p/2)))
      x11 <- y11 <- matrix(0, nobs, p)
      for(k in 1:p){
        screen(k)
        y11[, k] <- mgwerdiag$Lmaxr[[k]]/sum(mgwerdiag$Lmaxr[[k]])
        x11[, k] <- 1:length(object$fitted[!wzero])
        plot(x11[, k], y11[, k], xlab = "Index", ylab = "|Lmax|", ylim = c(0,1), ...) ; abline(2/nobs,0)
      }
      xx <- apply(x11, 2, list)#list(x5, x6, x7, x8)
      yy <- apply(y11, 2, list)#list(y5, y6, y7, y8)
    })
    if ((!(pick == 1))) {
      if (is.null(labels)) 
        labels <- names(object$lm$y)
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
    mgwerdiag
  else invisible()
}
