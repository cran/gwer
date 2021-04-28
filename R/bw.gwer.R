#' @title Optimization of Bandwidth for Geographically Weighted Elliptical Regression Model
#' @import spgwr
#' @import spdep
#' @importFrom GWmodel gw.dist gw.weight
#' @description The function compute the optimal bandwidth for a given geographically weighted elliptical regression using three differents methods: cross-validation, AIC and spatial validation. This optimal bandwidth optimzing the selected function.
#' @param formula regression model formula of a formula \code{object}.
#' @param family a description of the error distribution to be used in the model (see \code{family.elliptical} for more details of family functions).
#' @param data a SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}.
#' @param approach specified by CV for cross-validation approach, by AIC for corrected Akaike information criterion approach or by MI for spatial-validation approach.
#' @param spdisp if TRUE, by default, the dispersion parameter vary geographically in estimation process.
#' @param dispersion an optional fixed value for dispersion parameter.
#' @param kernel function chosen as follows:
#' gaussian: wgt = exp(-.5*(vdist/bw)^2);
#' exponential: wgt = exp(-vdist/bw);
#' bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
#' tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
#' boxcar: wgt=1 if dist < bw, wgt=0 otherwise.
#' @param adaptive if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance).
#' @param p	 the power of the Minkowski distance, default is 2 (Euclidean distance).
#' @param theta an angle in radians to rotate the coordinate system, default is 0
#' @param longlat if TRUE, great circle distances will be calculated.
#' @param dMat a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}.
#' @return returns the bandwidth optimization value.
#' @references Brunsdon, C., Fotheringham, A. S. and Charlton, M. E. (1996). 
#' Geographically weighted regression: a method for exploring spatial nonstationarity.
#' Geographical analysis, 28(4), 281-298. \doi{10.1111/j.1538-4632.1996.tb00936.x}
#' @references Cysneiros, F. J. A., Paula, G. A., and Galea, M. (2007). Heteroscedastic 
#' symmetrical linear models. Statistics & probability letters, 77(11), 1084-1090. 
#' \doi{10.1016/j.spl.2007.01.012} 
#' @references Fang, K. T., Kotz, S. and NG, K. W. (1990, ISBN:9781315897943).
#' Symmetric Multivariate and Related Distributions. London: Chapman and Hall.
#' @seealso \code{\link{gwer}}, \code{\link{elliptical}}, \code{\link{family.elliptical}}
#' @keywords Geographically weighted regression
#' @keywords Elliptical regression models
#' @keywords Bandwidth optimization
#' @examples
#' \donttest{
#' data(georgia, package="spgwr")
#' fit.formula <- PctBach ~ TotPop90 + PctRural + PctFB + PctPov
#' gwer.bw.n <- bw.gwer(fit.formula, data = gSRDF, family = Student(3), 
#'                     longlat = TRUE, adapt = TRUE)
#' }
#' @export

bw.gwer <- function(formula, family = Normal(), data, approach="CV", kernel="bisquare", 
                    adaptive = F, spdisp = 'local', dispersion, p=2, theta=0, longlat=F, 
                    dMat)
{
  if (is(data, "Spatial"))
  {
    dp.locat<-coordinates(data)
    data <- as(data, "data.frame")
  }
  else
  {
    stop("Given regression data must be Spatial*DataFrame")
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  dp.n <- nrow(data)
  
  if (missing(dMat))
  {
    dMat <- NULL
    DM.given<-F
    if(dp.n + dp.n <= 10000)
    {
      dMat <- gw.dist(dp.locat=dp.locat, rp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
      DM.given<-T
    }
  }
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }

  if(adaptive) {
    upper <- dp.n
    lower <- 20
  }
  else {
    if(DM.given)
    {
      upper<-range(dMat)[2]
      lower<-upper/5000
    }

    else
    {
      dMat<-NULL
      if (p==2)
      {
        b.box<-bbox(dp.locat)
        upper<-sqrt((b.box[1,2]-b.box[1,1])^2+(b.box[2,2]-b.box[2,1])^2)
        lower<-upper/5000
      }
      else
      {
        upper<-0
        for (i in 1:dp.n)
        {
          dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
          upper<-max(upper, range(dist.vi)[2])
        }
        lower<-upper/5000
      }
    }
  }
  
  if(spdisp == 'global' && missing(dispersion))
    stop ("Dispersion value not informed")
  
  if(spdisp == 'global')
    dispersion = rep(dispersion, dp.n)
  
  if(spdisp == 'local')
    dispersion = rep(NULL, dp.n)
  
  bw <- NA
  if (approach=="cv"||approach=="CV")
    bw <- gold2(gwer.cv, lower, upper, adapt.bw = adaptive, x, y, kernel, adaptive, dispersion, dp.locat, p, theta, longlat, dMat, family)
  else if (approach=="aic"||approach=="AIC")
    bw <- gold2(gwer.aic, lower, upper, adapt.bw = adaptive, x, y, kernel, adaptive, dispersion, dp.locat, p, theta, longlat, dMat, family)
  else if (approach=="aicc"||approach=="AICc")
    bw <- gold2(gwer.aicc, lower, upper, adapt.bw = adaptive, x, y, kernel, adaptive, dispersion, dp.locat, p, theta, longlat, dMat, family)
  else if (approach=="bic"||approach=="BIC")
    bw <- gold2(gwer.bic, lower, upper, adapt.bw = adaptive, x, y, kernel, adaptive, dispersion, dp.locat, p, theta, longlat, dMat, family)
  else if (approach=="mi"||approach=="MI")
    bw <- gold2(gwer.mi, lower, upper, adapt.bw = adaptive, x, y, kernel, adaptive, dispersion, dp.locat, p, theta, longlat, dMat, family)
  bw
}




gold2 <- function(fun, xL, xU, adapt.bw=F, ...){
  eps=1e-4 # working value - can be changed
  R <- (sqrt(5)-1)/2  # R <- 0.61803398....
  iter <- 1
  d <- R*(xU-xL)
  if (adapt.bw)
  {
    x1 <- floor(xL+d)
    x2 <- round(xU-d)
  }
  else
  {
    x1 <- xL+d
    x2 <- xU-d
  }  
  f1 <- eval(fun(x1,...))
  f2 <- eval(fun(x2,...))
  d1<-f2-f1
  # Establish initial value of xopt:
  if (f1 < f2)
    xopt  <-  x1
  else xopt  <-  x2
  # start main loop
  ea <- 100
  while ((abs(d) > eps) && (abs(d1) > eps)) {
    d <- R*d
    if   (f1 < f2) {
      xL  <-  x2
      x2  <-  x1
      if (adapt.bw)         
        x1 <- round(xL+d)
      else
        x1 <- xL+d
      #x1  <-  xL + d
      f2  <-  f1
      f1 <- eval(fun(x1,...))
    }
    else {
      xU  <-  x1
      x1  <-  x2
      if (adapt.bw)         
        x2 <- floor(xU - d)
      else
        x2  <-  xU - d 
      f1  <-  f2
      f2 <- eval(fun(x2,...))
    }
    iter  <-  iter + 1
    # Establish value of xopt after iteration:
    if    (f1 < f2)
      xopt  <-  x1
    else xopt  <-  x2
    d1<-f2-f1
  }
  xopt
}



gwer.cv <- function(bw, X, Y, kernel="bisquare", adaptive = F, dispersion, dp.locat, p = 2, theta = 0, longlat = F, dMat, family, verbose = T)
{
  dp.n <- length(dp.locat[,1])
  var.n <- ncol(X)
  if (is.null(dMat))
    DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }

  cv <- numeric(dp.n)

  for (i in 1:dp.n) {
    XX <- X[i, ]
    if (DM.given)
      dist.vi<-dMat[,i]
    else
    {
      dist.vi <- gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
    }

    w.i <- gw.weight(dist.vi, bw, kernel, adaptive)
    w.i[i] <- 0
    
    lm.i <- try(gwer.fit(Y = Y, X = X, gweights = w.i, family=family, offset = NULL, dispersion = dispersion[i], 
                         epsilon = 1e-04, maxit = 100, trace = F))
    if(!inherits(lm.i, "try-error") && lm.i$convergence == T) {
      #b <- coefficients(lm.i)
      yhat.noi<-XX%*%coefficients(lm.i)
      #cv[i] <- y[i] - (t(b) %*% xx)
      cv[i] <- Y[i]-yhat.noi
    }
    else{
      cv[i] <- Inf
      break
    }
  }
  if(!any(is.infinite(cv))){
    score <- t(cv) %*% cv
  }
  else {
    score <- Inf
  }
  
  if(verbose){
    if(adaptive)
      cat("Adaptive bandwidth (number of nearest neighbours)::", bw, "CV score:", score, "\n")
    else
      cat("Fixed bandwidth:", bw, "CV score:", score, "\n")
  }
#  if(is.nan(score))
#    score <- Inf
  score
}

gwer.aic <- function(bw, X, Y, kernel = "bisquare", adaptive = F, dispersion, dp.locat, p = 2, theta = 0, longlat = FALSE, dMat, family, verbose = T)
{
  dp.n<-length(dp.locat[,1])
  var.n <- ncol(X)
  if (is.null(dMat))
    DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }

  
  S <- matrix(nrow=dp.n, ncol=dp.n) ; flag <- 0 
  fittedgwr <- numeric(dp.n) ; dispersiongwr <- numeric(dp.n) 

  for (i in 1:dp.n) {
    if (DM.given)
      dist.vi<-dMat[,i]
    else
    {
      dist.vi <- gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
    }

    w.i <- gw.weight(dist.vi, bw, kernel, adaptive)

    lm.i <- try(gwer.fit(Y = Y, X = X, gweights = w.i, family=family, offset = NULL, dispersion = dispersion[i], 
                         epsilon = 1e-04, maxit = 100, trace = F))
    if(!inherits(lm.i, "try-error") && lm.i$convergence == T) {
      S[i,] <- lm.i$H[i,]
      fittedgwr[i] <- fitted.values(lm.i)[i]
      dispersiongwr[i] <- lm.i$dispersion
    } else {
      flag <- 1
      break
    }
  }
  
  if (flag == 0) {
    fittedgwr <- fittedgwr
    nu1 <- sum(diag(S)) ; res <- (Y - fittedgwr)/sqrt(dispersiongwr)
    logLik <- -0.5 * sum(log(dispersiongwr)) + sum(family$g0(res, df = family$df, s = family$s, 
                                                             r = family$r, alpha = family$alpha, mp = family$mp, 
                                                             epsi = family$epsi, sigmap = family$sigmap, k = family$k))
    score <- 2*nu1 - 2*logLik
  } else {
    score <- Inf
  }

  if(verbose){
    if(adaptive)
      cat("Adaptive bandwidth (number of nearest neighbours):", bw, "AIC value:", score, "\n")
    else
      cat("Fixed bandwidth:", bw, "AIC value:", score, "\n")
  }
#  if(is.nan(score))
#    score <- Inf
  score
} 

gwer.aicc <- function(bw, X, Y, kernel = "bisquare", adaptive = F, dispersion, dp.locat, p = 2, theta = 0, longlat = FALSE, dMat, family, verbose = T)
{
  dp.n<-length(dp.locat[,1])
  var.n <- ncol(X)
  if (is.null(dMat))
    DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }
  
  
  S <- matrix(nrow=dp.n, ncol=dp.n) ; flag <- 0 
  fittedgwr <- numeric(dp.n) ; dispersiongwr <- numeric(dp.n) 
  
  for (i in 1:dp.n) {
    if (DM.given)
      dist.vi<-dMat[,i]
    else
    {
      dist.vi <- gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
    }

    w.i <- gw.weight(dist.vi, bw, kernel, adaptive)
    
    lm.i <- try(gwer.fit(Y = Y, X = X, gweights = w.i, family=family, offset = NULL, dispersion = dispersion[i], 
                          epsilon = 1e-04, maxit = 100, trace = F))
    if(!inherits(lm.i, "try-error") && lm.i$convergence == T) {
      S[i,] <- lm.i$H[i,]
      fittedgwr[i] <- fitted.values(lm.i)[i]
      dispersiongwr[i] <- lm.i$dispersion
    } else {
      flag <- 1
      break
    }
  }
  
  if (flag == 0) {
    fittedgwr <- fittedgwr
    nu1 <- sum(diag(S)) ; res <- (Y - fittedgwr)/sqrt(dispersiongwr)
    logLik <- -0.5 * sum(log(dispersiongwr)) + sum(family$g0(res, df = family$df, s = family$s, 
                                                             r = family$r, alpha = family$alpha, mp = family$mp, 
                                                             epsi = family$epsi, sigmap = family$sigmap, k = family$k))
    score <- (dp.n*(dp.n + nu1))/(dp.n-2-nu1) - 2*logLik 	#2*nu1 + 2*(nu1)*(nu1+1)/(dp.n-nu1-1) - 2*logLik
  } else {
    score <- Inf
  }
  
  if(verbose){
    if(adaptive)
      cat("Adaptive bandwidth (number of nearest neighbours):", bw, "AICc value:", score, "\n")
    else
      cat("Fixed bandwidth:", bw, "AICc value:", score, "\n")
  }
  #  if(is.nan(score))
  #    score <- Inf
  score
}

gwer.bic <- function(bw, X, Y, kernel = "bisquare", adaptive = F, dispersion, dp.locat, p = 2, theta = 0, longlat = FALSE, dMat, family, verbose = T)
{
  dp.n<-length(dp.locat[,1])
  var.n <- ncol(X)
  if (is.null(dMat))
    DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }
  
  
  S <- matrix(nrow=dp.n, ncol=dp.n) ; flag <- 0 
  fittedgwr <- numeric(dp.n) ; dispersiongwr <- numeric(dp.n) 
  
  for (i in 1:dp.n) {
    if (DM.given)
      dist.vi<-dMat[,i]
    else
    {
      dist.vi <- gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
    }

    w.i <- gw.weight(dist.vi, bw, kernel, adaptive)
    
    lm.i <- try(gwer.fit(Y = Y, X = X, gweights = w.i, family=family, offset = NULL, dispersion = dispersion[i], 
                          epsilon = 1e-04, maxit = 100, trace = F))
    if(!inherits(lm.i, "try-error") && lm.i$convergence == T) {
      S[i,] <- lm.i$H[i,]
      fittedgwr[i] <- fitted.values(lm.i)[i]
      dispersiongwr[i] <- lm.i$dispersion
    } else {
      flag <- 1
      break
    }
  }
  
  if (flag == 0) {
    fittedgwr <- fittedgwr
    nu1 <- sum(diag(S)) ; res <- (Y - fittedgwr)/sqrt(dispersiongwr)
    logLik <- -0.5 * sum(log(dispersiongwr)) + sum(family$g0(res, df = family$df, s = family$s, 
                                                             r = family$r, alpha = family$alpha, mp = family$mp, 
                                                             epsi = family$epsi, sigmap = family$sigmap, k = family$k))
    score <- log(dp.n)*nu1 - 2*logLik
  } else {
    score <- Inf
  }
  
  if(verbose){
    if(adaptive)
      cat("Adaptive bandwidth (number of nearest neighbours):", bw, "BIC value:", score, "\n")
    else
      cat("Fixed bandwidth:", bw, "BIC value:", score, "\n")
  }
  #  if(is.nan(score))
  #    score <- Inf
  score
} 


gwer.mi <- function(bw, X, Y, kernel="bisquare", adaptive = F, dispersion, dp.locat, p=2, theta = 0, longlat=FALSE, dMat, family, verbose = T)
{
  dp.n<-length(dp.locat[,1])
  var.n <- ncol(X)
  if (is.null(dMat))
    DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
      stop ("Dimensions of dMat are not correct")
  }

  residuals <- numeric(dp.n) ; resid <- numeric(dp.n) ; flag <- 0 
  s <- numeric(dp.n) ; rs  <- numeric(dp.n) ; S <- matrix(nrow=dp.n, ncol=dp.n)

  for (i in 1:dp.n) {
    if (DM.given)
      dist.vi<-dMat[,i]
    else
    {
      dist.vi <- gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
    }

    w.i <- gw.weight(dist.vi, bw, kernel, adaptive)

    lm.i <- try(gwer.fit(Y = Y, X = X, gweights = w.i, family=family, offset = NULL, dispersion = dispersion[i], 
                         epsilon = 1e-04, maxit = 100, trace = F))
    if(!inherits(lm.i, "try-error") && lm.i$convergence == T) {
      Xd <- lm.i$Xmodel
      resid <- lm.i$residuals
      H <- lm.i$Hat[i,]
      h <- diag(H)/(lm.i$scalevariance * lm.i$scale)
      rs <- resid/sqrt(lm.i$scalevariance * (1 - h))
      residuals[i] <- rs[i]
    } else {
      flag <- 1
      break
    }
  }
  
  if (flag == 0) {
    distcoord <- knn2nb(knearneigh(dp.locat, longlat=longlat))
    col.test <- nb2listw(distcoord, style="W")
    morani <- moran.test(residuals, col.test, alternative = 'two.sided')
    score  <- abs(as.numeric(morani$estimate[1]))
  } else {
    score <- Inf
  }
  
  if(verbose){
    if(adaptive)
      cat("Adaptive bandwidth (number of nearest neighbours):", bw, "Moran I.:", score, "\n")
    else
      cat("Fixed bandwidth:", bw, "Moran I.:", score, "\n")
  }
  if(is.nan(score))
    score <- Inf
  score
}