#' @title Geographically Weighted Elliptical Regression
#' @import spgwr
#' @import graphics
#' @import methods
#' @import sp
#' @import spData
#' @import maptools
#' @importFrom GWmodel gw.dist gw.weight
#' @name gwer
#' @aliases print.gwer
#' @description The function fit geographically weighted elliptical regression  model to explore the non-stationarity for a certain bandwidth and  weighting function.
#' @param formula regression model formula as in \code{glm}.
#' @param bandwidth value of the selected bandwidth used in the weighting function (see \code{bw.gwer} for bandwidth optimization).
#' @param kernel function chosen as follows:
#' gaussian: wgt = exp(-.5*(vdist/bw)^2);
#' exponential: wgt = exp(-vdist/bw);
#' bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
#' tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
#' boxcar: wgt=1 if dist < bw, wgt=0 otherwise
#' @param p the power of the Minkowski distance, default is 2, i.e. the Euclidean distance
#' @param theta an angle in radians to rotate the coordinate system, default is 0
#' @param dMat a pre-specified distance matrix, it can be calculated by the function gw.dist
#' @param regression.points a Spatial*DataFrame object, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp; Note that no diagnostic information will returned if it is assigned.
#' @param adapt defines the type of bandwidth used. either NULL (default) or a proportion between 0 and 1 of observations to include in weighting scheme (k-nearest neighbours).
#' @param hatmatrix if TRUE, return the hatmatrix as a component of the result.
#' @param spdisp if TRUE dispersion parameter varies geographically.
#' @param family a description of the error distribution to be used in the model (see \code{\link{family.elliptical}} for details of family functions).
#' @param data model data frame, or may be a SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}.
#' @param dispersion an optional fixed value for dispersion parameter.
#' @param weights an optional numeric vector of weights to be used in the fitting process.
#' @param subset an optional numeric vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain NAs (see \code{glm}).
#' @param method the method to be used in fitting local models. The default method "bw.gwer" uses Fisher's scoring method. The alternative "model.frame" returns the model frame and does no fitting.
#' @param longlat TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers. If x is a SpatialPoints object, the value is taken from the object itself.
#' @param control a list of parameters for controlling the fitting process. For \code{elliptical} this is passed by \code{glm.control}.
#' @param model a logical value indicating whether model frame should be included as a component of the return.
#' @param x a logical value indicating whether the response vector used in the fitting process should be returned as components of the return.
#' @param y a logical value indicating whether model matrix used in the fitting process should be returned as components of the return.
#' @param contrasts an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.
#' @param parplot if TRUE the parameters boxplots are plotted.
#' @param offset this can be used to specify an a priori known component to be included in the linear predictor during fitting as in \code{glm}.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return returns an object of class \dQuote{gwer}, a list with follow components:   
#' \item{SDF}{a SpatialPointsDataFrame (may be gridded) or SpatialPolygonsDataFrame object (see package \pkg{sp}) with fit.points, weights, GWR coefficient estimates, dispersion and the residuals in its \code{data} slot.}
#' \item{coef}{the matrices of coefficients, standard errors and significance values for parameters hypothesis test.}
#' \item{dispersion}{either the supplied argument or the estimated dispersion with standard error.}
#' \item{hat}{hat matrix of the geographically weighted elliptical model.}
#' \item{lm}{elliptical global regression on the same model formula.}  
#' \item{results}{a list of results values for fitted geographically weighted elliptical model.}  
#' \item{bandwidth}{the bandwidth used in geographical weighting function.}
#' \item{fitted}{the fitted mean values of the geographically weighted elliptical model.}
#' \item{hatmatrix}{a logical value indicating if hatmatrix was considered}
#' \item{gweights}{a matrix with the geographical weighting for all local elliptical models.}
#' \item{family}{the \code{family} object used.}
#' \item{flm}{a matrix with the fitted values for all local elliptical models.}
#' \item{adapt}{the \code{adapt} object used.}
#' \item{kernel}{the \code{kernel} object used.}
#' \item{spdisp}{the \code{spdisp} object used.}
#' \item{this.call}{the function call used.}
#' \item{longlat}{the \code{longlat} object used.}
#' @references Brunsdon, C., Fotheringham, A. S. and Charlton, M. E. (1996). 
#' Geographically weighted regression: a method for exploring spatial nonstationarity.
#' Geographical analysis, 28(4), 281-298. \doi{10.1111/j.1538-4632.1996.tb00936.x}
#' @references Cysneiros, F. J. A., Paula, G. A., and Galea, M. (2007). Heteroscedastic 
#' symmetrical linear models. Statistics & probability letters, 77(11), 1084-1090. 
#' \doi{10.1016/j.spl.2007.01.012} 
#' @references Fang, K. T., Kotz, S. and NG, K. W. (1990, ISBN:9781315897943).
#' Symmetric Multivariate and Related Distributions. London: Chapman and Hall.
#' @seealso \code{\link{bw.gwer}}, \code{\link{elliptical}}, \code{\link{family.elliptical}}
#' @keywords Geographically weighted regression
#' @keywords Bandwidth optimization
#' @keywords Elliptical model
#' @examples
#' \donttest{
#' data(georgia, package = "spgwr")
#' fit.formula <- PctBach ~ TotPop90 + PctRural + PctFB + PctPov
#' gwer.bw.t <- bw.gwer(fit.formula, data = gSRDF, family = Student(3), adapt = TRUE)
#' gwer.fit.t <- gwer(fit.formula, data = gSRDF, family = Student(3), bandwidth = gwer.bw.t, 
#'                    adapt = TRUE, parplot = FALSE, hatmatrix = TRUE, spdisp = TRUE, 
#'                    method = "gwer.fit")
#' print(gwer.fit.t) 
#' }
#' @rdname gwer
#' @export

gwer <- function (formula, data, regression.points, bandwidth, kernel = "bisquare", 
                  p = 2, theta = 0, adapt = NULL, hatmatrix = FALSE, family = Normal, longlat = NULL, dMat,
                  weights, dispersion = NULL, subset, na.action = "na.fail", method = "gwer.fit",
                  control = glm.control(epsilon = 1e-04, maxit = 100, trace = F), model = FALSE,
                  x = FALSE, y = TRUE, contrasts = NULL, offset, spdisp = TRUE, parplot = FALSE,...)
{
  call <- match.call()
  this.call <- match.call()
  dist <- as.character(call$family)[1]
  user.def <- F
  
  if (!charmatch(dist, c("Normal", "Cauchy", "Student", "Gstudent", "LogisI", "LogisII", "Glogis", "Cnormal", "Powerexp"),nomatch = F)) 
    dist <- match.arg(dist, c("Normal", "Cauchy", "Student", "Gstudent", "LogisI", "LogisII", "Glogis", "Cnormal", "Powerexp"))
  else user.def <- T
  
  resid_name <- "residuals"
  p4s <- as.character(NA)
  Polys <- NULL
  if (is(data, "SpatialPolygonsDataFrame")) 
    Polys <- as(data, "SpatialPolygons")
  
  if (missing(regression.points)) {
    rp.given <- FALSE
    regression.points <- data
    rp.locat <- coordinates(data)
    hatmatrix <- T
  }
  else {
    rp.given <- TRUE
    hatmatrix <- F
    if (is(regression.points, "Spatial")) {
      rp.locat <- coordinates(regression.points)
    }
    else if (is.numeric(regression.points) && dim(regression.points)[2] == 
             2) 
      rp.locat <- regression.points
    else {
      warning("Output loactions are not packed in a Spatial object,and it has to be a two-column numeric vector")
      rp.locat <- dp.locat
    }
  }
  if (is(data, "Spatial")) {
    p4s <- proj4string(data)
    dp.locat <- coordinates(data)
    data <- as(data, "data.frame")
  }
  else {
    stop("Given regression data must be Spatial*DataFrame")
  }
  
  if (is.null(longlat) || !is.logical(longlat)) 
    longlat <- FALSE
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data)) 
    data <- environment(formula)

  
  if (!charmatch(method, c("model.frame", "gwer.fit"), F)) 
    stop(paste("\n unimplemented method:", method))
  m <- match.call(expand.dots = FALSE)
  m$family <- m$method <- m$control <- m$model <- m$dispersion <- m$x <- m$y <- m$bandwidth  <- m$kernel <- m$adapt <- m$longlat <-m$contrasts <- m$offset <- m$hatmatrix <- m$spdisp <- m$parplot <- m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  dp.n <- length(model.extract(m, "response"))
  if (method == "model.frame") 
    return(m)
  {
    if (!missing(family) && !charmatch(dist, c("Normal", "Cauchy", "Student", 
                                               "Gstudent", "LogisI", "LogisII", 
                                               "Glogis", "Cnormal", "Powerexp"), F)) 
      cat(paste("\n work with user-defined family:", call$family, "\n"))
  }
  if (!missing(dispersion) && is.numeric(dispersion) && !(dispersion > 0)) 
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
  
  
  method <- "gwer.fit"
  elliptical.fitter <- get(method)
  offset4fit <- offset
  
  fit <- elliptical.fitter(X = X, Y = Y, offset = offset4fit, family = family, 
                           dispersion = dispersion, maxit = control$maxit, 
                           epsilon = control$epsilon, trace = control$trace, ...)
  if (!spdisp)
    dispersiongwr <- dispersion <- fit$dispersion
  if (hatmatrix) 
    se.fit <- TRUE
  
  n <- NROW(rp.locat)
  rownames(rp.locat) <- NULL
  if (is.null(colnames(rp.locat))) 
    colnames(rp.locat) <- c("x", "y")
  p <- NCOL(X)

  if (NROW(X) != NROW(regression.points)) 
    stop("Input data and coordinates have different dimensions")

  DM.given <- F
  if (missing(dMat)) {
    DM.given <- F
    DM1.given <- F
    if (dp.n + n <= 5000) {
      dMat <- gw.dist(dp.locat = dp.locat, rp.locat = rp.locat, 
                      p = p, theta = theta, longlat = longlat)
      DM.given <- T
    }
    else {
      dMat <- matrix(0, 1, 1)
    }
  }
  else {
    DM.given <- T
    DM1.given <- T
    dim.dMat <- dim(dMat)
    if (dim.dMat[1] != dp.n || dim.dMat[2] != n) 
      stop("Dimensions of dMat are not correct")
  }

  v_resids <- numeric(n) ; v_fitteds <- numeric(n)

  if (any(wzero)) {
    nas <- is.na(fit$coef)
    fitted[wpos] <- fit$fitted.values/w[wpos]
    fitted[wzero] <- X.org[wzero, !nas] %*% as.vector(fit$coefficients[!nas]) + 
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
    nas <- is.na(fit$coefficients)
    pasgn <- asgn[!nas]
    if (df.residuals > 0) 
      fit$assign.residual <- (rank + 1):length(Y)
    fit$R.assign <- pasgn
    fit$x.assign <- asgn
  }
  
  gweights <- matrix(0,n,n) ; local.fitted <- matrix(0,n,n) ; sum.w <- numeric(n) ; resid.ord <- matrix(0,n,1)
    

  if(spdisp){ 
    dispersiongwr <- numeric(n) 
    gwr.b <- matrix(nrow = nobs, ncol = ncol(X)+1)
    stderror <- matrix(0,n,p+1) ; zvalue <- matrix(0,n,p) ; pvalue <- matrix(0,n,p)
    colnames(gwr.b) <- c(colnames(X), "dispersion")
    colnames(stderror) <- c(paste(colnames(X),rep(".se",p), sep = ""), "dispersion.se")
  } else {
    gwr.b <- matrix(nrow = nobs, ncol = ncol(X))
    stderror <- matrix(0,n,p) ; zvalue <- matrix(0,n,p) ; pvalue <- matrix(0,n,p)
    colnames(gwr.b) <- c(colnames(X))    
    colnames(stderror) <- paste(colnames(X),rep(".se",p), sep = "")
  }
  colnames(zvalue) <- colnames(pvalue) <- c(colnames(X)) 
  
  if (hatmatrix)
    Hat <- matrix(0,n,n)
  
  for (i in 1:n) {
    if (DM.given) 
      dist.vi <- dMat[, i]
    else {
      if (rp.given) 
        dist.vi <- gw.dist(dp.locat, rp.locat, focus = i, 
                           p, theta, longlat)
      else dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, 
                              p = p, theta = theta, longlat = longlat)
    }

    w.i <- gw.weight(dist.vi, bandwidth, kernel, adapt)
    if (any(w.i < 0 | is.na(w.i))) 
      stop(paste("Invalid weights for i:", i))
    lm.i <- elliptical.fitter(X = X, Y = Y, gweights = w.i, offset = offset4fit, 
                              family = family, dispersion = dispersion, 
                              maxit = control$maxit, epsilon = control$epsilon, 
                              trace = control$trace, ...)
    sum.w[i] <- sum(w.i)
    v_fitteds[i] <- lm.i$fitted.values[i]
    if(spdisp){ 
      dispersiongwr[i] <- lm.i$dispersion
      gwr.b[i, ] <- c(coefficients(lm.i), dispersiongwr[i])
    } else{
      gwr.b[i, ] <- coefficients(lm.i)
    }
    
    if (any(wzero)) {
      nas <- is.na(lm.i$coefficients)
      fitted.i[wpos] <- lm.i$fitted.values/w[wpos]
      fitted.i[wzero] <- X.org[wzero, !nas] %*% as.vector(lm.i$coefficients[!nas]) + 
        if (length(offset.org) > 1) 
          offset.org[wzero]
      else 0
      lm.i$fitted.values <- fitted.i
      resid.i[wpos] <- lm.i$resid
      resid.i[wzero] <- (Y.org[wzero] - fitted.i[wzero])/sqrt(lm.i$dispersion)
      lm.i$residuals <- resid.i
      q1.i[wpos] <- lm.i$q1
      q2.i[wpos] <- lm.i$q2
      q1.i[wzero] <- family$g1(resid[wzero], df = family$df, 
                               alpha = family$alpha, mp = family$mp, epsi = family$epsi, 
                               sigmap = family$sigmap, k = family$k)
      q2.i[wzero] <- -2 * q1.i[wzero]
      lm.i$q1 <- q1
      lm.i$q2 <- q2
    }
    else lm.i$fitted.values <- lm.i$fitted.values/w
    lm.i$weights <- w
    lm.i$gweights <- w.i
    gweights[i,] <- lm.i$gweights
    names(lm.i$fitted.values) <- names(lm.i$residuals) <- names(lm.i$q1) <- names(lm.i$q2) <- NULL
    rank.i <- lm.i$rank
    df.r <- length(if (exists("X.org", frame = sys.nframe())) Y.org else Y) - 
      rank.i - sum(w == 0)
    asgn.i <- attr(if (exists("X.org", frame = sys.nframe())) X.org else X, 
                   "assign")
    if (rank < p) {
      nas.i <- is.na(lm.i$coefficients)
      pasgn.i <- asgn[!nas]
      if (df.residuals > 0) 
        lm.i$assign.residual <- (rank.i + 1):length(Y)
      lm.i$R.assign <- pasgn.i
      lm.i$x.assign <- asgn.i
    }
    lm.i <- c(lm.i, list(assign = asgn.i, df.residuals = df.residuals, 
              family = family, user.def = user.def, formula = as.vector(attr(Terms, 
              "formula")), terms = Terms, contrasts = attr(X, 
              "contrasts"), control = control, call = call))
    
    if (y) 
      lm.i$y <- if (exists("Y.org", frame = sys.nframe())) 
        Y.org
    else Y
    names(lm.i$y) <- NULL
    if (x) 
      lm.i$X <- if (exists("X.org", frame = sys.nframe())) 
        X.org
    else X
    if (model) 
      lm.i$model <- m
    attr(lm.i, "class") <- c("elliptical", "glm", "lm")
    
    local.fitted[i,] <- lm.i$fitted.values

    R <- lm.i$R[(1:p), (1:p)]
    covun <- solve(qr(R))
    rowlen <- sqrt(diag(covun))

    if (hatmatrix) {
      Hat[i,] <- lm.i$Hat[i,] ; nu1 <- sum(diag(Hat)) ; nu2 <- sum(diag(t(Hat)%*%Hat)) ; alpha.pv <- (nu1/p)
      if(spdisp){
        stderror[i, ] <- c(rowlen[1:p] * sqrt(lm.i$dispersion/lm.i$scale), sqrt(4*lm.i$dispersion^2/(sum(lm.i$gweights)*lm.i$scaledispersion)))
        zvalue[i, ] <- c(gwr.b[i, 1:p])/stderror[i, 1:p]
        pvalue[i, ] <- 2 * pnorm(-abs(zvalue[i, ]))#/(nu1/p)
      } else {
        stderror[i, ] <- rowlen[1:p] * sqrt(lm.i$dispersion/lm.i$scale)
        zvalue[i, ] <- c(gwr.b[i, ])/stderror[i, ]
        pvalue[i, ] <- 2 * pnorm(-abs(zvalue[i, ]))#/(nu1/p)
      }
    }
  }
  
  if (any(wzero)) {
    fittedgwr[wpos] <- v_fitteds/w[wpos]
    fittedgwr[wzero] <- X.org[wzero, !nas] %*% as.vector(fit$coefficients[!nas]) + 
      if (length(offset.org) > 1) 
        offset.org[wzero]
    else 0
    residualgwr[wpos] <- resid.ord
    residualgwr[wzero] <- (Y.org[wzero] - fitted[wzero])/sqrt(dispersiongwr)
  }	
  else {
    fittedgwr <- v_fitteds/w
    residualgwr <- resid.ord
  }
  
  df <- data.frame(sum.w = sum.w, gwr.b)
  if (hatmatrix)
    df <- data.frame(sum.w = sum.w, gwr.b, stderror)
  df[[resid_name]] <- v_resids
  SDF <- SpatialPointsDataFrame(coords = rp.locat, data = df, 
                                proj4string = CRS(p4s))

  if (hatmatrix) {
    res <- (Y - fittedgwr)/sqrt(dispersiongwr)
    logLik <- -0.5 * sum(log(dispersiongwr)) + sum(family$g0(res, df = family$df, s = family$s, r = family$r,
                                                             alpha = family$alpha, mp = family$mp, epsi = family$epsi,
                                                             sigmap = family$sigmap, k = family$k))
    AIC <- 2*nu1 - 2*logLik
    AICc <- 2*nu1 - 2*logLik + 2*(nu1)*(nu1+1)/(dp.n-nu1-1)
    BIC <- log(dp.n)*nu1 - 2*logLik
    
    edf <- dp.n - 2 * nu1 + nu2
    
    results <- list(df.residuals = df.residuals, nu1 = nu1, nu2 = nu2, edf = edf, 
                    logLik = - 2*logLik, AIC = AIC, AICc = AICc, BIC = BIC)
    
    local.coef <- list(est = gwr.b, se = stderror, pvalue = pvalue)
  } else {
    local.coef <- list(est = gwr.b)
    results <- Hat <- NULL
  }
  
  if(parplot == TRUE){
    par(mfrow = c(ceiling(sqrt(p)),ceiling(sqrt(p))))
    for(i in 1:p){
      boxplot(gwr.b[ , i], ylab=substitute(expression(beta[j]), list(j=i-1)))
    }
    if(spdisp)
      boxplot(dispersiongwr, ylab=expression(phi))
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

  z <- list(SDF = SDF, coef = local.coef, dispersion = dispersiongwr, hat = Hat, lm = fit, 
            results = results, bandwidth = bandwidth, fitted = fittedgwr, residual = residualgwr, 
            hatmatrix = hatmatrix, gweights = gweights, family = fit$family, flm = local.fitted, 
            adapt = adapt, kernel = kernel, spdisp = spdisp, alpha.pv = alpha.pv, this.call = this.call, 
            longlat = longlat, control = control, offset = offset)
  
  class(z) <- "gwer"
  invisible(z)
}




gwer.fit <- function (X, Y, gweights=NULL, offset, family, dispersion, 
                      maxit, epsilon, trace) 
{
  n <- nrow(X)
  if(is.null(gweights))
    gweights <- rep(1,n)
  if (is.null(offset)) 
    offset <- rep(0, n)
  
  p <- ncol(X)
  aux.model <- glm.fit(x = X, y = Y, offset = offset, weights = gweights,
                       family = gaussian())
  attr(aux.model, "class") <- c("glm", "lm")
  start <- aux.model$coef
  
  is.null.disp <- is.null(dispersion)
  elliptical.disp <- !is.null.disp && !is.numeric(dispersion)
  if (is.null.disp){
    options(warn = -1)
    dispersion <- (summary(aux.model)$dispersion)
    options(warn = 0) 
  }
  if (elliptical.disp){
    options(warn = -1)
    dispersion <- (summary(aux.model)$dispersion)
    options(warn = 0)
  }
  args <- resid(aux.model)/sqrt(dispersion)
  
  if (any(nas <- is.na(start))) {
    names(nas) <- dimnames(X)[[2]]
    X <- X[, !nas]
    aux.model <- glm.fit(x = X, y = Y, offset = offset, weights = gweights,
                         family = gaussian())
    attr(aux.model, "class") <- c("glm", "lm")
    start <- aux.model$coef
    options(warn = -1)
    dispersion <- (summary(aux.model)$dispersion)
    options(warn = 0)
  }
  
  linearconv <- TRUE
  iter <- 1
  error1 <- error2 <- error3 <- 0
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
      w.h <- as.vector(-2 * w.1 * gweights)
      
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
      new.dispersion <- sum((-2 * w.1 * gweights) * abs.res^2)/(sum(gweights))
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
  if (maxit > 1 && iter == maxit){
    linearconv <- F
    warning(paste("\n linear convergence not obtained in", 
                  maxit, "iterations"))
  }
  coefs <- rep(NA, length(nas))
  coefs[!nas] <- start
  names(coefs) <- names(nas)
  names(dispersion) <- "dispersion"
  
  fitted <- as.vector(X %*% start + offset)
  
  residuals <- (Y - fitted)/sqrt(dispersion)
  w.1 <- family$g1(residuals, df = family$df, s = family$s, 
                   r = family$r, alpha = family$alpha, mp = family$mp, 
                   epsi = family$epsi, sigmap = family$sigmap, k = family$k)
  w.2 <- -2 * w.1 * gweights
  if (any(w.2 < 0)) 
    cat("\n --- negative iterative weights returned! --- \n")
  
  if (is.null.disp) {
    Xd <- cbind(X, residuals) ; rank <- dim(X)[2]
    Rnames <- dimnames(X)[[2]]
    dimnames(Xd)[[2]] <- c(Rnames, "scale")
  } else {
    Xd <- cbind(X) ; rank <- dim(X)[2]
    Rnames <- dimnames(X)[[2]]
    dimnames(Xd)[[2]] <- Rnames
  }

  nn <- is.null(Rnames)
  Rnames <- list(dimnames(Xd)[[2]], dimnames(Xd)[[2]])
  R <- t(Xd)%*%diag(as.vector(gweights))%*%Xd
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
  w.2.h <- as.vector(-2 * w.1 * gweights)
  Hat <- X %*% solve(t(X) %*% diag(w.2.h) %*% X) %*% t(X) %*% diag(w.2.h)
  C <- solve(t(X) %*% diag(w.2.h) %*% X) %*% t(X) %*% diag(w.2.h)
  H <- X %*% solve(t(X) %*% diag(as.vector(gweights)) %*% X) %*% t(X) %*% diag(as.vector(gweights))
  Haux <- solve(t(X) %*% diag(as.vector(gweights)) %*% X) %*% t(X) %*% diag(as.vector(gweights))
  
  fit <- list(coefficients = coefs, dispersion = dispersion, gweights = gweights,
              fixed = !is.null.disp, residuals = residuals, fitted.values = fitted, 
              loglik = loglik, convergence = linearconv, Hat = Hat, Ci = C, H = H, Haux = Haux,
              Wg = family$g1(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), Wgder = family$g5(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), v = -2 * family$g1(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), rank = rank, R = as.matrix(R),  iter = iter - 
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
              "Powerexp", F)) family$k, Xmodel = matrix(Xd[, (1:rank)],nrow(Xd), rank))
  fit
}


#' @rdname gwer
#' @method print gwer
#' @noRd
#' @export
print.gwer <- function(x, ...) {
  if(class(x) != "gwer") 
    stop("not a gwer object")
  cat("Call:\n")
  print(x$this.call)
  cat("Kernel function:", x$kernel, "\n")
  n <- length(x$lm$residuals)
  if (is.null(x$adapt)) cat("Fixed bandwidth:", x$bandwidth, "\n")
  else cat("Adaptive quantile: ", x$adapt, " (about ", floor(x$adapt*n), " of ", n, " data points)\n", sep="")
  m <- length(x$lm$coefficients)

  if(!x$spdisp){
    df0 <- as(x$SDF, "data.frame")[,(1+(1:m)), drop=FALSE]
  } else {
    df0 <- as(x$SDF, "data.frame")[,(1+(1:(m+1))), drop=FALSE]
  }
  if (any(is.na(df0))) {
    df0 <- na.omit(df0)
    warning("NAs in coefficients dropped")
  }
  CM <- t(apply(df0, 2, summary))[,c(1:3,5,6)]
  if (is.null(dim(CM))) CM <- t(as.matrix(CM))
  if(!x$spdisp){
     CM <- cbind(CM, c(coefficients(x$lm)))
  } else {
    CM <- cbind(CM, c(coefficients(x$lm),x$lm$dispersion))
  }
  colnames(CM) <- c(colnames(CM)[1:5], "Global")
  printCoefmat(CM)
  cat("Number of data points:", n, "\n")
  if(!x$spdisp)
    cat("Global dispersion:", x$lm$dispersion, "\n")
  if (x$hatmatrix) {
    cat("Effective number of parameters (residual: 2traceS - traceS'S):", 2*x$results$nu1 -
          x$results$nu2, "\n")
    cat("Effective degrees of freedom (residual: 2traceS - traceS'S):", x$results$edf, "\n")
    cat("Effective number of parameters (model: traceS):",
        x$results$nu1, "\n")
    cat("Effective degrees of freedom (model: traceS):",
        (n - x$results$nu1), "\n")
    cat("AIC:", x$results$AIC, " AICc:", x$results$AICc, " BIC:", x$results$BIC, "\n")
    #cat("Residual sum of squares:", x$results$rss, "\n")
    #cat("Quasi-global R2:", (1 - (x$results$rss/x$gTSS)), "\n")
  }
  invisible(x)
}
