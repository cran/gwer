#' @title Optimization of Bandwidth for Geographically Weighted Elliptical Regression Model
#' @import stats
#' @import spdep
#' @description The function compute the optimal bandwidth for a given geographically weighted elliptical regression using three differents methods: cross-validation, AIC and spatial validation. This optimal bandwidth optimzing the selected function.
#' @param formula regression model formula as in \code{glm}.
#' @param data model data frame, or may be a SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}.
#' @param coords matrix of coordinates of points representing the spatial positions of the observations.
#' @param adapt defines the type of bandwidth used. Either TRUE: find the proportion between 0 and 1 of observations to include in weighting scheme (k-nearest neighbours) or FALSE: find global bandwidth.
#' @param gweight geographical weighting function, at present \code{gwr.Gauss()} default.
#' @param method  type of the method used to the compute of residuals. Is \code{cv} for drop-1 cross-validation (default), \code{aic} for AIC optimisation (depends on assumptions about AIC degrees of freedom) or \code{sv} for spatial validation.
#' @param verbose if TRUE (default) reports the progress of search for bandwidth.
#' @param longlat TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if x is a SpatialPoints object, the value is taken from the object itself.
#' @param family a description of the error distribution to be used in the model (see \code{family.elliptical} for more details of family functions).
#' @param RMSE default FALSE to correspond with CV scores in newer references (sum of squared CV errors), if TRUE the previous behaviour of scoring by LOO CV RMSE.
#' @param weights an optional numeric vector of weights to be used in the fitting process, beware of scaling issues. Only used with the cross-validation method, probably unsafe.
#' @param tol the desired accuracy to be passed to \code{optimize}.
#' @param show.error.messages default FALSE. may be set to TRUE to see error messages if \code{gwer.sel} returns without a value.
#' @param maxit maximum number of iterations in model fit
#' @return returns the bandwidth optimization value.
#' @references Brunsdon, C., Fotheringham, A. S. and Charlton, M. E. (1996). 
#' Geographically weighted regression: a method for exploring spatial nonstationarity.
#' Geographical analysis, 28(4), 281-298. \url{https://doi.org/10.1111/j.1538-4632.1996.tb00936.x}
#' @references Cysneiros, F. J. A., Paula, G. A., and Galea, M. (2007). Heteroscedastic 
#' symmetrical linear models. Statistics & probability letters, 77(11), 1084-1090. 
#' \url{https://doi.org/10.1016/j.spl.2007.01.012} 
#' @references Fang, K. T., Kotz, S. and NG, K. W. (1990, ISBN:9781315897943).
#' Symmetric Multivariate and Related Distributions. London: Chapman and Hall.
#' @seealso \code{\link{gwer}}, \code{\link{elliptical}}, \code{\link{family.elliptical}}
#' @keywords Geographically weighted regression
#' @keywords Bandwidth optimization
#' @keywords Elliptical model
#' @examples
#' data(columbus, package="spData")
#' gwer.bw <- gwer.sel(CRIME ~ INC, data=columbus, family = Normal(),
#'                  coords=cbind(columbus$X, columbus$Y))
#' \donttest{
#' data(columbus, package="spData")
#' gwer.bw <- gwer.sel(CRIME ~ INC, data=columbus, family = Student(df=4),
#'                  coords=cbind(columbus$X, columbus$Y), method = "aic")
#' }
#' @export

gwer.sel <- function(formula, data = list(), coords, adapt=FALSE, gweight=gwr.Gauss, weights,
                    method="cv", verbose=TRUE, longlat=NULL, family = Normal() , RMSE=FALSE,  
                    tol=.Machine$double.eps^0.25, show.error.messages=FALSE, maxit = 100) 
{
  if (!is.logical(adapt)) 
    stop("adapt must be logical")
  if (is(data, "Spatial")) {
    if (!missing(coords)) 
      warning("data is Spatial* object, ignoring coords argument")
    coords <- coordinates(data)
    if (is.null(longlat) || !is.logical(longlat)) {
      if (!is.na(is.projected(data)) && !is.projected(data)) {
        longlat <- TRUE
      }
      else {
        longlat <- FALSE
      }
    }
    data <- as(data, "data.frame")
  }
  
  if (is.null(longlat) || !is.logical(longlat)) 
    longlat <- FALSE
  if (missing(coords)) 
    stop("Observation coordinates have to be given")
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  dp.n <- length(model.extract(mf, "response"))
  weights <- as.vector(model.extract(mf, "weights"))
  # set up default weights
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  if (is.null(weights)) weights <- rep(as.numeric(1), dp.n)
  if (any(is.na(weights))) stop("NAs in weights")
  if (any(weights < 0)) stop("negative weights")
  
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  
  if(method == "cv"){
    bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
    difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
    if (any(!is.finite(difmin))) 
      difmin[which(!is.finite(difmin))] <- 0
    
    if (!adapt){
      beta1 <- difmin/1000
      beta2 <- difmin
      opt <- optimize(gwer.cv.f, lower = beta1, upper = beta2, y = y, x = x, 
                      maximum = FALSE, family = family, coords = coords, gweight = gweight, 
                      verbose = verbose, longlat = longlat, RMSE = RMSE, weights=weights, 
                      tol = tol, show.error.messages=show.error.messages, maxit = maxit)
    }
    else{
      beta1 <- 0
      beta2 <- 1
      opt <- optimize(gwer.cv.adapt.f, lower = beta1, upper = beta2, y = y , x = x,
                      maximum = FALSE, family = family, coords = coords, gweight = gweight,
                      verbose = verbose, longlat = longlat, RMSE = RMSE, weights=weights, 
                      tol = tol, show.error.messages=show.error.messages, maxit = maxit)
    }
    bdwt <- opt$minimum
    res <- bdwt
  }
  
  if(method == "aic"){	
    bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
    difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
    if (any(!is.finite(difmin))) 
      difmin[which(!is.finite(difmin))] <- 0

    if (!adapt){
      beta1 <- difmin/1000
      beta2 <- difmin
      opt <- optimize(gwer.aic.f, lower = beta1, upper = beta2, y = y, x = x, 
                      maximum = FALSE, family = family, coords = coords, gweight = gweight, 
                      verbose = verbose, longlat = longlat, RMSE = RMSE, weights=weights, 
                      tol = tol, show.error.messages=show.error.messages, maxit = maxit)
    }
    else{
      beta1 <- 0
      beta2 <- 1
      opt <- optimize(gwer.aic.adapt.f, lower = beta1, upper = beta2, y = y , x = x,
                      maximum = FALSE, family = family, coords = coords, gweight = gweight,
                      verbose = verbose, longlat = longlat, RMSE = RMSE, weights=weights, 
                      tol = tol, show.error.messages=show.error.messages, maxit = maxit)
    }
    bdwt <- opt$minimum
    res <- bdwt
  }
  
  if(method == "mi"){	
    bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
    difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
    if (any(!is.finite(difmin))) 
      difmin[which(!is.finite(difmin))] <- 0

    if (!adapt){
      beta1 <- difmin/1000
      beta2 <- difmin
      opt <- optimize(gwer.mi.f, lower = beta1, upper = beta2, y = y, x = x, 
                      maximum = FALSE, family = family, coords = coords, gweight = gweight, 
                      verbose = verbose, longlat = longlat, RMSE = RMSE, weights=weights, 
                      tol = tol, show.error.messages=show.error.messages, maxit = maxit)
    }
    else{
      beta1 <- 0
      beta2 <- 1
      opt <- optimize(gwer.mi.adapt.f, lower = beta1, upper = beta2, y = y , x = x,
                      maximum = FALSE, family = family, coords = coords, gweight = gweight,
                      verbose = verbose, longlat = longlat, RMSE = RMSE, weights=weights, 
                      tol = tol, show.error.messages=show.error.messages, maxit = maxit)
    }
    bdwt <- opt$minimum
    res <- bdwt
  }
  
  res
}



gwer.cv.f <- function(bandwidth, y, x, coords, gweight, family, verbose=TRUE, longlat=FALSE, 
                      RMSE=FALSE, weights = NULL, show.error.messages=FALSE, maxit = 100)
{
  n <- NROW(x)
  #    	m <- NCOL(x)
  if (is.null(weights)) 
    weights <- rep(1, n)
  cv <- numeric(n)
  options(show.error.messages = show.error.messages)
  for (i in 1:n) {
    xx <- x[i, ]
    dxs <- spDistsN1(coords, coords[i,], longlat=longlat)
    if (!is.finite(dxs[i])) dxs[i] <- .Machine$double.xmax/2
    w.i <- gweight(dxs^2, bandwidth)
    #		w.i <- gweight(spDistsN1(coords, coords[i,], longlat=longlat)^2, bandwidth)
    w.i[i] <- 0
    w.i <- w.i * weights
    if (any(w.i < 0 | is.na(w.i)))
      stop(paste("Invalid weights for i:", i))
    lm.i <- try(gwer.fit(Y = y, X = x, gweights = w.i, family=family, offset = NULL, dispersion = NULL, maxit = maxit, 
                         epsilon = 1e-04, trace = F))
    if(!inherits(lm.i, "try-error") && lm.i$convergence == T) {
      b <- coefficients(lm.i)
      cv[i] <- y[i] - (t(b) %*% xx)
    }
    else{
      cv[i] <- Inf
      break
    }
  }
  if (!any(is.infinite(cv))){
    score <- sum(t(cv) %*% cv)
    if (RMSE) score <- sqrt(score/n)
    #    	score <- sqrt(sum(t(cv) %*% cv)/n)
  }
  else {
    score <- Inf
  }

  if (!show.error.messages) options(show.error.messages = TRUE)
  if (verbose) cat("Fixed bandwidth:", bandwidth, "CV score:", score, "\n")
  score
}


gwer.cv.adapt.f <- function(q, y, x, coords, gweight, family, verbose=TRUE, longlat=FALSE, 
                            RMSE=FALSE, weights = NULL, show.error.messages=FALSE, maxit = 100)
{
  n <- NROW(x)
  #    	m <- NCOL(x)
  if (is.null(weights)) 
    weights <- rep(1, n)
  cv <- numeric(n)
  bw <- gw.adapt(dp=coords, fp=coords, quant=q, longlat=longlat)
  options(show.error.messages = show.error.messages)
  for (i in 1:n) {
    xx <- x[i, ]
    dxs <- spDistsN1(coords, coords[i,], longlat=longlat)
    if (!is.finite(dxs[i])) dxs[i] <- .Machine$double.xmax/2
    w.i <- gweight(dxs^2, bw[i])
    w.i[i] <- 0
    w.i <- w.i * weights
    if (any(w.i < 0 | is.na(w.i)))
      stop(paste("Invalid weights for i:", i))
    lm.i <- try(gwer.fit(Y = y, X = x, gweights = w.i, family=family, offset = NULL, dispersion = NULL, maxit = maxit, 
                         epsilon = 1e-04, trace = F))
    if(!inherits(lm.i, "try-error") && lm.i$convergence == T) {
      b <- coefficients(lm.i)
      cv[i] <- y[i] - (t(b) %*% xx)
    }
    else{
      cv[i] <- Inf
      break
    }
  }

  if (!any(is.infinite(cv))){
    score <- sum(t(cv) %*% cv)
    if (RMSE) score <- sqrt(score/n)
    #    	score <- sqrt(sum(t(cv) %*% cv)/n)
  }
  else {
    score <- Inf
  }
  
  if (!show.error.messages) options(show.error.messages = TRUE)
  if (verbose) cat("Adaptive bandwidth:", q, "CV score:", score, "\n")
  score
}



gwer.aic.f <- function(bandwidth, y, x, coords, gweight, family, verbose=TRUE, longlat=FALSE, 
                       RMSE=FALSE, weights = NULL, show.error.messages=FALSE, maxit = 100)
{
  n <- NROW(x)
  #    	m <- NCOL(x)
  if (is.null(weights)) 
    weights <- rep(1, n)
  Hat <- matrix(nrow=n, ncol=n) ; flag <- 0 
  fittedgwr <- numeric(n) ; dispersiongwr <- numeric(n) 
  options(show.error.messages = show.error.messages)
  for (i in 1:n) {
    #        	xx <- x[i, ]
    dxs <- spDistsN1(coords, coords[i,], longlat=longlat)
    if (!is.finite(dxs[i])) dxs[i] <- .Machine$double.xmax/2
    w.i <- gweight(dxs^2, bandwidth)
    #		w.i <- gweight(spDistsN1(coords, coords[i,], longlat=longlat)^2, bandwidth)
    if (any(w.i < 0 | is.na(w.i)))
      stop(paste("Invalid weights for i:", i))
    lm.i <- try(gwer.fit(Y = y, X = x, gweights = w.i, family=family, offset = NULL, dispersion = NULL, maxit = maxit, 
                         epsilon = 1e-04, trace = F))
    if(!inherits(lm.i, "try-error") && lm.i$convergence == T) {
      fittedgwr[i] <- fitted.values(lm.i)[i]
      Hat[i,] <- t(x[i,]) %*% solve(t(x) %*% diag(w.i) %*% x) %*% t(x) %*% diag(w.i)
      dispersiongwr[i] <- lm.i$dispersion
    } else {
      flag <- 1
    }
  }
  
  if (flag == 0) {
    fittedgwr <- fittedgwr/weights
    nu1 <- sum(diag(Hat)) ; res <- (y - fittedgwr)/sqrt(dispersiongwr)
    logLik <- -0.5 * sum(log(dispersiongwr)) + sum(family$g0(res, df = family$df, s = family$s, 
                                                             r = family$r, alpha = family$alpha, mp = family$mp, 
                                                             epsi = family$epsi, sigmap = family$sigmap, k = family$k))
    score <- 2*nu1 - 2*logLik + 2*(nu1)*(nu1+1)/(n-nu1-1)	#-2*logLik + (n*(n + nu1)/(n - 2 - nu1))
  } else {
    score <- Inf
  }
  if (!show.error.messages) options(show.error.messages = TRUE)
  if (verbose) cat("Fixed bandwidth:", bandwidth, "AICc:", score, "\n")
  score
} 



gwer.aic.adapt.f <- function(q, y, x, coords, gweight, family, verbose=TRUE, longlat=FALSE, 
                             RMSE=FALSE, weights = NULL, show.error.messages=FALSE, maxit = 100)
{
  n <- NROW(x)
  #    	m <- NCOL(x)
  bw <- gw.adapt(dp=coords, fp=coords, quant=q, longlat=longlat)
  Hat <- matrix(nrow=n, ncol=n) ; flag <- 0 
  fittedgwr <- numeric(n) ; dispersiongwr <- numeric(n) 
  options(show.error.messages = show.error.messages)
  for (i in 1:n) {
    #       	xx <- x[i, ]
    dxs <- spDistsN1(coords, coords[i,], longlat=longlat)
    if (!is.finite(dxs[i])) dxs[i] <- .Machine$double.xmax/2
    w.i <- gweight(dxs^2, bw[i])
    #		w.i <- gweight(spDistsN1(coords, coords[i,], longlat=longlat)^2, bw[i])
    if (any(w.i < 0 | is.na(w.i)))
      stop(paste("Invalid weights for i:", i))
    lm.i <- try(gwer.fit(Y = y, X = x, gweights = w.i, family=family, offset = NULL, dispersion = NULL, maxit = maxit, 
                         epsilon = 1e-04, trace = F))
    if(!inherits(lm.i, "try-error") && lm.i$convergence == T) {
      fittedgwr[i] <- fitted.values(lm.i)[i]
      Hat[i,] <- t(x[i,]) %*% solve(t(x) %*% diag(w.i) %*% x) %*% t(x) %*% diag(w.i)
      dispersiongwr[i] <- lm.i$dispersion
    } else {
      flag <- 1
    }
  }
  
  if (flag == 0) {
    fittedgwr <- fittedgwr/weights
    nu1 <- sum(diag(Hat)) ; res <- (y - fittedgwr)/sqrt(dispersiongwr)
    logLik <- -0.5 * sum(log(dispersiongwr)) + sum(family$g0(res, df = family$df, s = family$s, 
                                                             r = family$r, alpha = family$alpha, mp = family$mp, 
                                                             epsi = family$epsi, sigmap = family$sigmap, k = family$k))
    score <- 2*nu1 - 2*logLik + 2*(nu1)*(nu1+1)/(n-nu1-1)	#-2*logLik + (n*(n + nu1)/(n - 2 - nu1))
  } else {
    score <- Inf
  }
  if (!show.error.messages) options(show.error.messages = TRUE)
  if (verbose) cat("Adaptive bandwidth:", q, "AICc:", score, "\n")
  score
}



gwer.mi.f <- function(bandwidth, y, x, coords, gweight, family, verbose=TRUE, longlat=FALSE, 
                      RMSE=FALSE, weights = NULL, show.error.messages=FALSE, maxit = 100)
{
  n <- NROW(x)
  #    	m <- NCOL(x)
  if (is.null(weights)) 
    weights <- rep(1, n)
  wzero <- (weights == 0)
  residuals <- numeric(n) ; resid <- numeric(n) ; flag <- 0 
  h <- numeric(n) ; rs  <- numeric(n) ; H <- matrix(nrow=n, ncol=n)
  options(show.error.messages = show.error.messages)
  for (i in 1:n) {
    #        	xx <- x[i, ]
    dxs <- spDistsN1(coords, coords[i,], longlat=longlat)
    if (!is.finite(dxs[i])) dxs[i] <- .Machine$double.xmax/2
    w.i <- gweight(dxs^2, bandwidth)
    #		w.i <- gweight(spDistsN1(coords, coords[i,], longlat=longlat)^2, bandwidth)
    if (any(w.i < 0 | is.na(w.i)))
      stop(paste("Invalid weights for i:", i))
    lm.i <- try(gwer.fit(Y = y, X = x, gweights = w.i, family=family, offset = NULL, dispersion = NULL, maxit = maxit, 
                         epsilon = 1e-04, trace = F))
    if(!inherits(lm.i, "try-error") && lm.i$convergence == T) {
      Xd <- diag(c(weights[!wzero])) %*% lm.i$Xmodel[!wzero, ]
      resid <- lm.i$residuals
      H <- Xd %*% solve(t(Xd) %*% diag(w.i) %*% Xd) %*% t(Xd) %*% diag(w.i)
      h <- diag(H)/(lm.i$scalevariance * lm.i$scale)
      rs <- resid/sqrt(lm.i$scalevariance * (1 - h))
      residuals[i] <- rs[i]
    } else {
      flag <- 1
    }
  }
  
  if (flag == 0) {
    distcoord <- knn2nb(knearneigh(coords, longlat=longlat))
    col.test <- nb2listw(distcoord, style="W")
    morani <- moran.test(residuals, col.test, alternative = 'two.sided')
    score  <- abs(as.numeric(morani$estimate[1]))
  } else {
    score <- Inf
  }
  
  if (!show.error.messages) options(show.error.messages = TRUE)
  if (verbose) cat("Fixed bandwidth:", bandwidth, "Moran I.:", score, "\n")
  score
}


gwer.mi.adapt.f <- function(q, y, x, coords, gweight, family, verbose=TRUE, longlat=FALSE, 
                            RMSE=FALSE, weights = NULL, show.error.messages=TRUE, maxit = 100)
{
  n <- NROW(x)
  #    	m <- NCOL(x)
  if (is.null(weights)) 
    weights <- rep(1, n)
  wzero <- (weights == 0)
  residuals <- numeric(n) ; resid <- numeric(n) ; flag <- 0 
  h <- numeric(n) ; rs  <- numeric(n) ; H <- matrix(nrow=n, ncol=n)
  bw <- gw.adapt(dp=coords, fp=coords, quant=q, longlat=longlat)
  options(show.error.messages = show.error.messages)
  for (i in 1:n) {
    #       	xx <- x[i, ]
    dxs <- spDistsN1(coords, coords[i,], longlat=longlat)
    if (!is.finite(dxs[i])) dxs[i] <- .Machine$double.xmax/2
    w.i <- gweight(dxs^2, bw[i])
    #		w.i <- gweight(spDistsN1(coords, coords[i,], longlat=longlat)^2, bw[i])
    if (any(w.i < 0 | is.na(w.i)))
      stop(paste("Invalid weights for i:", i))
    lm.i <- try(gwer.fit(Y = y, X = x, gweights = w.i, family=family, offset = NULL, dispersion = NULL, maxit = maxit, 
                         epsilon = 1e-04, trace = F))
    if(!inherits(lm.i, "try-error") && lm.i$convergence == T) {
      Xd <- diag(c(weights[!wzero])) %*% lm.i$Xmodel[!wzero, ]
      resid <- lm.i$residuals
      H <- Xd %*% solve(t(Xd) %*% diag(w.i) %*% Xd) %*% t(Xd) %*% diag(w.i)
      h <- diag(H)/(lm.i$scalevariance * lm.i$scale)
      rs <- resid/sqrt(lm.i$scalevariance * (1 - h))
      residuals[i] <- rs[i]
    } else {
      flag <- 1
    }
  }
  
  if (flag == 0) {
    distcoord <- knn2nb(knearneigh(coords, longlat=longlat))
    col.test <- nb2listw(distcoord, style="W")
    morani <- moran.test(residuals, col.test, alternative = 'two.sided')
    score  <- abs(as.numeric(morani$estimate[1]))
  } else {
    score <- Inf
  }
  
  if (!show.error.messages) options(show.error.messages = TRUE)
  if (verbose) cat("Adaptive bandwidth:", q, "Moran I.:", score, "\n")
  score
}
