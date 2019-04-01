#' @title Geographically Weighted Elliptical Regression
#' @import spgwr
#' @import stats
#' @import graphics
#' @import methods
#' @import sp
#' @import spData
#' @import maptools
#' @description The function implements geographically weighted elliptical regression to explore the non-stationarity for certain global bandwidth and chosen weighting scheme.
#' @param formula regression model formula as in \code{glm}.
#' @param coords matrix of coordinates of points representing the spatial positions of the observations.
#' @param bandwidth value of the selected bandwidth (see \code{gwer.sel} for bandwidth optimization).
#' @param gweight geographical weighting function, at present \code{gwr.Gauss()} is default.
#' @param adapt defines the type of bandwidth used. either NULL (default) or a proportion between 0 and 1 of observations to include in weighting scheme (k-nearest neighbours).
#' @param spdisp if TRUE dispersion parameter varies geographically.
#' @param family a description of the error distribution to be used in the model (see \code{elliptical.family} for details of family functions).
#' @param data model data frame, or may be a SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}.
#' @param dispersion an optional fixed value for dispersion parameter.
#' @param weights an optional vector of weights to be used in the fitting process for local models.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param fit.points an object containing the coordinates of fit points, often an object from package \pkg{sp}. If missing, the coordinates given through the data argument object, or the coords argument are used.
#' @param na.action a function which indicates what should happen when the data contain NAs (see \code{glm}).
#' @param method the method to be used in fitting local models. The default method "gwer.fit" uses Fisher's scoring method. The alternative "model.frame" returns the model frame and does no fitting.
#' @param longlat TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers. If x is a SpatialPoints object, the value is taken from the object itself.
#' @param control a list of parameters for controlling the fitting process. For \code{gwer.fit} this is passed to \code{glm.control}.
#' @param model a logical value indicating whether model frame should be included as a component of the returned value.
#' @param x a logical value indicating whether the response vector used in the fitting process should be returned as components of the returned value.
#' @param y a logical value indicating whether model matrix used in the fitting process should be returned as components of the returned value.
#' @param contrasts an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.
#' @param parplot if TRUE the parameters boxplots are plotted.
#' @param offset this can be used to specify an a priori known component to be included in the linear predictor during fitting as in \code{glm}.
#' @param type character that indicates the type of residuals should consider as return.
#' @param gwr.diag if TRUE is calculated the diagnostic measures of the model and provided in return.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return A list of class \dQuote{gwer}:   
#' \item{SDF}{a SpatialPointsDataFrame (may be gridded) or SpatialPolygonsDataFrame object (see package \pkg{sp}) with fit.points, weights, GWR coefficient estimates, dispersion and the residuals of \code{type} in its \code{data} slot.}
#' \item{coef}{regression parameters matrix of the fitted model.}
#' \item{se}{standard errors matrix for the parameters of the fitted model.}
#' \item{pvalue}{p-value matrix for the significance tests of parameters of the fitted model.}
#' \item{lhat}{hat matrix of the geographically weighted elliptical model.}
#' \item{lm}{elliptical regression on the same model formula.}  
#' \item{Weights}{matrix for geographical weighting.}
#' \item{results}{a list of results values for fitted geographically weighted elliptical model.}
#' \item{fitted}{the fitted mean values of the geographically weighted elliptical model.}
#' \item{diag}{a list of diagnostic matrices ('leverage', 'global influence' and 'local influence'), see \code{elliptical.diag} for more details.}
#' \item{residuals}{a list of all residuals type ('ordinal', 'studentized' and 'deviance').}
#' \item{this.call}{the function call used.}
#' @references Brunsdon, C., Fotheringham, A. S. and Charlton, M. E. (1996). 
#' Geographically weighted regression: a method for exploring spatial nonstationarity.
#' Geographical analysis, 28(4), 281-298. \url{https://doi.org/10.1111/j.1538-4632.1996.tb00936.x}
#' @references Fang, K. T., Kotz, S. and NG, K. W. (1990, ISBN:9781315897943).
#' Symmetric Multivariate and Related Distributions. London: Chapman and Hall.
#' @seealso \code{\link{gwer.sel}}, \code{\link{elliptical}}, \code{\link{family.elliptical}}
#' @keywords spatial
#' @keywords elliptical model
#' @examples
#' data(columbus, package="spData")
#' fit.lm <- lm(CRIME ~ INC, data=columbus)
#' summary(fit.lm)
#' gwer.bw <- gwer.sel(CRIME ~ INC, data=columbus, family = Normal(),
#'                  coords=cbind(columbus$X, columbus$Y))
#' fit.gwer <- gwer(CRIME ~ INC, family = Normal(), bandwidth = gwer.bw, 
#'                  parplot = TRUE, data=columbus, method = "gwer.fit",
#'                  coords=cbind(columbus$X, columbus$Y))
#' \donttest{
#' data(columbus, package="spData")
#' fit.lm <- lm(CRIME ~ INC, data=columbus)
#' summary(fit.lm)
#' gwer.bw <- gwer.sel(CRIME ~ INC, data=columbus, family = Student(df=4),
#'                  coords=cbind(columbus$X, columbus$Y), method = 'aic')
#' gwer.fitt <- gwer(CRIME ~ INC, family = Student(df=4), bandwidth = gwer.bw, 
#'                  parplot = TRUE, data=columbus, method = "gwer.fit",
#'                  coords=cbind(columbus$X, columbus$Y))
#' }
#' @export

gwer <- function (formula, coords, bandwidth, gweight = gwr.Gauss, adapt = NULL, spdisp = F,
          family = Normal, data = sys.parent(), dispersion = NULL, weights, subset,
          fit.points, na.action = "na.fail", method = "gwer.fit", longlat = NULL,
          control = glm.control(epsilon = 1e-04, maxit = 100, trace = F), model = F,
          x = F, y = T, contrasts = NULL, parplot = F, offset, type='pearson', gwr.diag = F,...) 
{
  call <- match.call()
  this.call <- match.call()
  dist <- as.character(call$family)[1]
  user.def <- F
  
  if (!charmatch(dist, c("Normal", "Cauchy", "Student", "Gstudent", 
                         "LogisI", "LogisII", "Glogis", "Cnormal", "Powerexp"), 
                 nomatch = F)) 
    dist <- match.arg(dist, c("Normal", "Cauchy", "Student", 
                              "Gstudent", "LogisI", "LogisII", "Glogis", "Cnormal", 
                              "Powerexp"))
  else user.def <- T
  
  resid_name <- paste(type, "resids", sep = "_")
  p4s <- as.character(NA)
  Polys <- NULL
  if (is(data, "SpatialPolygonsDataFrame")) 
    Polys <- as(data, "SpatialPolygons")
  
  if (is(data, "Spatial")) {
    if (!missing(coords)) 
      warning("data is Spatial* object, ignoring coords argument")
    coords <- coordinates(data)
    p4s <- proj4string(data)
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
  if (is.null(colnames(coords))) 
    colnames(coords) <- c("coord.x", "coord.y")
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
  #if (!is.null(gwr.diag) && gwr.diag != 'idex' && gwr.diag != 'mean') 
  #  stop("diagnostic method not implemented (choose between 'idex' or 'mean')")
  
  
  if (!charmatch(method, c("model.frame", "gwer.fit"), F)) 
    stop(paste("\n unimplemented method:", method))
  m <- match.call(expand.dots = F)
  m$family <- m$method <- m$control <- m$model <- m$dispersion <- m$x <- m$y <- m$bandwidth  <- m$gweight <- m$adapt <- m$longlat <-m$contrasts <- m$offset <- m$gwr.diag <- m$parplot <- m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  if (method == "model.frame") 
    return(m)
  {
    if (!missing(family) && !charmatch(dist, c("Normal", 
                                               "Cauchy", "Student", "Gstudent", "LogisI", "LogisII", 
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
  
  fit <- elliptical.fitter(X = X, Y = Y,  offset = offset4fit, family = family, 
                           dispersion = dispersion, maxit = control$maxit, 
                           epsilon = control$epsilon, trace = control$trace, ...)
  
  if (missing(fit.points)) {
    fp.given <- FALSE
    fit.points <- coords
    colnames(fit.points) <- colnames(coords)
  }
  else fp.given <- TRUE
  griddedObj <- FALSE
  if (is(fit.points, "Spatial")) {
    Polys <- NULL
    if (is(fit.points, "SpatialPolygonsDataFrame")) {
      Polys <- Polygons(fit.points)
      fit.points <- coordinates(fit.points)
    }
    else {
      griddedObj <- gridded(fit.points)
      fit.points <- coordinates(fit.points)
    }
  }
  n <- NROW(fit.points)
  rownames(fit.points) <- NULL
  if (is.null(colnames(fit.points))) 
    colnames(fit.points) <- c("x", "y")
  p <- NCOL(x)

  if (NROW(X) != NROW(coords)) 
    stop("Input data and coordinates have different dimensions")
  
  if (is.null(adapt)) {  
    if (!missing(bandwidth)) {
      bw <- bandwidth
      bandwidth <- rep(bandwidth, nobs)
    }
    else stop("Bandwidth must be given for non-adaptive weights")
  }
  else {
    bandwidth <- gw.adapt(dp = coords, fp = fit.points, quant = adapt, 
                          longlat = longlat)
    bw <- bandwidth
  }
  
  if (any(bandwidth < 0)) 
    stop("Invalid bandwidth")
  gwr.b <- matrix(nrow = nobs, ncol = ncol(X))
  v_resids <- numeric(n)
  v_fitteds <- numeric(n)
  colnames(gwr.b) <- colnames(X)
  
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
  
  
  lhat <- NA ; gWeights <- matrix(0,n,n) ; Vi <- matrix(0,n,n)
  stderror <- matrix(0,n,p+1) ; zvalue <- matrix(0,n,p+1) ; pvalue <- matrix(0,n,p+1)
  colnames(stderror) <- c(paste(colnames(X),rep(".se",p), sep = ""), "dispersion.se") ; colnames(zvalue) <- colnames(pvalue) <- c(colnames(X),"dispersion")
  sum.w <- numeric(n) ; dispersiongwr <- numeric(n) ; Hat <- matrix(0,n,n)
  residgwr <- numeric(n) ; resid.diagst <- resid.diagcd <- matrix(0,n,1)
  if(!is.null(gwr.diag)){
    GL <- GLbeta <- GLphi <- Lmaxr <- Cic <- Cih <- DGbeta <- DGphi <- hat <- matrix(0,n,n)
    Cmaxc <- list() ; for(j in 1:p){Cmaxc[[j]] <- matrix(0,n,n)}
  }
  
  for (i in 1:n) {
    dxs <- spDistsN1(coords, fit.points[i, ], longlat = longlat)
    if (any(!is.finite(dxs))) 
      dxs[which(!is.finite(dxs))] <- .Machine$double.xmax/2
    w.i <- gweight(dxs^2, bandwidth[i])
    if (any(w.i < 0 | is.na(w.i))) 
      stop(paste("Invalid weights for i:", i))
    lm.i <- elliptical.fitter(X = X, Y = Y, gweights = w.i, offset = offset4fit, 
                              family = family, dispersion = dispersion, 
                              maxit = control$maxit, epsilon = control$epsilon, 
                              trace = control$trace, ...)
    sum.w[i] <- sum(w.i)
    gwr.b[i, ] <- coefficients(lm.i)
    if (!fp.given) 
      v_resids[i] <- residuals.elliptical(lm.i, type = type)[i]
    else is.na(v_resids[i]) <- TRUE
    
    v_fitteds[i] <- fitted.values(lm.i)[i]
    
    if (any(wzero)) {
      nas <- is.na(fit$coef)
      fitted.i[wpos] <- lm.i$fitted.values/w[wpos]
      fitted.i[wzero] <- X.org[wzero, !nas] %*% as.vector(fit$coef[!nas]) + 
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
    gWeights[i,] <- lm.i$gweights
    names(lm.i$fitted.values) <- names(lm.i$residuals) <- names(lm.i$q1) <- names(lm.i$q2) <- NULL
    rank.i <- lm.i$rank
    df.r <- length(if (exists("X.org", frame = sys.nframe())) Y.org else Y) - 
      rank.i - sum(w == 0)
    asgn.i <- attr(if (exists("X.org", frame = sys.nframe())) X.org else X, 
                   "assign")
    if (rank < p) {
      nas.i <- is.na(lm.i$coef)
      pasgn.i <- asgn[!nas]
      if (df.residuals > 0) 
        lm.i$assign.residual <- (rank.i + 1):length(Y)
      lm.i$R.assign <- pasgn.i
      lm.i$x.assign <- asgn.i
    }
    lm.i <- c(lm.i, list(assign = asgn, df.residuals = df.residuals, 
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
    
    
    if(spdisp == TRUE)
    	dispersiongwr[i] <- lm.i$dispersion
    else
    	dispersiongwr[i] <- fit$dispersion
    
    
    R <- lm.i$R[(1:p), (1:p)]
    covun <- solve(qr(R))
    rowlen <- sqrt(diag(covun))
    
    stderror[i, ] <- c(rowlen[1:p] * sqrt(lm.i$dispersion/lm.i$scale), sqrt(4*lm.i$dispersion^2/(sum(lm.i$gweights)*lm.i$scaledispersion)))
    zvalue[i, ] <- c(gwr.b[i, ], dispersiongwr[i])/stderror[i, ]
    pvalue[i, ] <- 2 * pnorm(-abs(zvalue[i, ]))
    
    
    Vi[i,] <- lm.i$v ; Hat[i,] <- X[i,] %*% solve(t(X) %*% diag(w.i) %*% X) %*% t(X) %*% diag(w.i)
    if(!is.null(gwr.diag)){
      diag.i <- gwer.diag(lm.i)
      
      resid.diagst[i] <- diag.i$rs[i] ; resid.diagcd[i] <- diag.i$rd[i]
      GL[i,] <- diag.i$GL ; GLbeta[i,] <- diag.i$GLbeta ; GLphi[i,] <- diag.i$GLphi
      Cih[i,] <- diag.i$Cih ; Cic[i,] <- diag.i$Cic ; Lmaxr[i,] <- diag.i$Lmaxr
      for(j in 1:p){Cmaxc[[j]][i,] <- diag.i$Cmaxc[j,]}
      DGbeta[i,] <- diag.i$DGbeta ; DGphi[i,] <- diag.i$DGphi ; hat[i,] <- diag.i$h
    }
  }
  
  if (any(wzero)) {
    fittedgwr[wpos] <- v_fitteds/w[wpos]
    fittedgwr[wzero] <- X.org[wzero, !nas] %*% as.vector(fit$coef[!nas]) + 
      if (length(offset.org) > 1) 
        offset.org[wzero]
    else 0
    residgwr[wpos] <- v_resids
    residgwr[wzero] <- (Y.org[wzero] - fitted[wzero])/sqrt(fit$dispersiongwr)
  }	
  else {
    fittedgwr <- v_fitteds/w
    residgwr <- v_resids
  }
  
  df <- data.frame(sum.w = sum.w, gwr.b, dispersion = dispersiongwr, stderror)
  df[[resid_name]] <- residgwr
  SDF <- SpatialPointsDataFrame(coords = fit.points, data = df, 
                                proj4string = CRS(p4s))
  if (griddedObj) {
    gridded(SDF) <- TRUE
  }
  else {
    if (!is.null(Polys)) {
      df <- data.frame(SDF@data)
      rownames(df) <- sapply(slot(Polys, "polygons"), function(i) slot(i, "ID"))
      SDF <- SpatialPolygonsDataFrame(Sr = Polys, data = df)
    }
  }
  
  
  
  fit <- c(fit, list(assign = asgn, df.residuals = df.residuals, 
                     family = family, user.def = user.def, formula = as.vector(attr(Terms, 
                                                                                    "formula")), terms = Terms, contrasts = attr(X, 
                                                                                                                                 "contrasts"), control = control, call = call))
  
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
  
  nu1 <- sum(diag(Hat)) ; res <- (Y - fittedgwr)/sqrt(dispersiongwr)
  logLik <- -0.5 * sum(log(dispersiongwr)) + sum(family$g0(res, df = family$df, s = family$s, 
                                                           r = family$r, alpha = family$alpha, mp = family$mp, 
                                                           epsi = family$epsi, sigmap = family$sigmap, k = family$k))
  
  
  AIC <- 2*nu1 - 2*logLik
  AICc <- 2*nu1 - 2*logLik + 2*(nu1)*(nu1+1)/(n-nu1-1)
  BIC <- log(n)*nu1 - 2*logLik
  
  results <- list(df.residuals = df.residuals, dispersion = dispersiongwr, V = Vi, logLik = - 2*logLik, 
                  AIC = AIC, AICc = AICc, BIC = BIC)
  
  diagnostic <- NULL
  if(gwr.diag == TRUE){
    residgwr <- list(ord = residgwr, st = resid.diagst, cd = resid.diagcd)
    diagnostic <- list(GL = diag(GL), GLbeta = diag(GLbeta), GLphi = diag(GLphi), 
                       Cih = diag(Cih), Cic = diag(Cic), Lmaxr = diag(Lmaxr),
                       Cmaxc = lapply(Cmaxc,function(x){diag(x)}), DGbeta = diag(DGbeta),
                       DGphi = diag(DGphi))
  }

  
  if(parplot == TRUE){
    par(mfrow = c(ceiling(sqrt(p)),ceiling(sqrt(p))))
    for(i in 1:p){
      boxplot(gwr.b[ , i], ylab=substitute(expression(beta[j]), list(j=i-1)))
    }
    boxplot(dispersiongwr, ylab=expression(phi))
  }
  
  z <- list(SDF = SDF, coef = gwr.b, se = stderror, pvalue = pvalue, lhat = Hat, lm = fit, Weights = gWeights, results = results, fitted = fittedgwr, diag = diagnostic,
            bandwidth = bw, adapt = adapt, hatmatrix = FALSE, gweight = deparse(substitute(gweight)), residuals = residgwr,
            fp.given = fp.given,  this.call = this.call)
  
  class(z) <- "gwr"
  invisible(z)
}




gwer.fit <- function (X, Y, gweights=NULL, offset, family, dispersion, 
                      maxit, epsilon, trace) 
{
  n <- nrow(X)
  if(is.null(gweights))
    gweights <- rep(1,n)
  if (is.null(offset)) {
    offset <- rep(0, n)
  }
  
  p <- ncol(X)
  aux.model <- glm.fit(x = X, y = Y, offset = offset,
                       family = gaussian())
  attr(aux.model, "class") <- c("glm", "lm")
  start <- aux.model$coef
  
  is.null.disp <- is.null(dispersion)
  elliptical.disp <- !is.null.disp && !is.numeric(dispersion)
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
  
  
  linearconv <- TRUE
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
      new.dispersion <- mean((-2 * w.1 * gweights) * abs.res^2)
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
    rank <- dim(X)[2]
    Rnames <- dimnames(X)[[2]]
    Xd <- cbind(X, residuals)
  }
  dimnames(Xd)[[2]] <- c(Rnames, "scale")
  nn <- is.null(Rnames)
  Rnames <- list(dimnames(Xd)[[2]], dimnames(Xd)[[2]])
  R <- t(Xd) %*% diag(gweights) %*% Xd
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
  fit <- list(coefficients = coefs, dispersion = dispersion, gweights = gweights,
              fixed = !is.null.disp, residuals = residuals, fitted.values = fitted, 
              loglik = loglik, convergence = linearconv, Wg = family$g1(residuals, df = family$df, 
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


gwer.diag <- function (ellipticalfit, weighting = "observed") 
{
  Wi <- diag(ellipticalfit$gweights)
  scalevariance <- ellipticalfit$scalevariance
  scaledispersion <- ellipticalfit$scaledispersion
  scale <- ellipticalfit$scale
  family <- ellipticalfit$family
  user.def <- ellipticalfit$user.def
  f.name <- family[[1]]
  dispersion <- ellipticalfit$dispersion
  w <- if (is.null(ellipticalfit$weights)) 
    rep(1, length(ellipticalfit$residuals))
  else ellipticalfit$weights
  wzero <- (w == 0)
  
  resid <- ellipticalfit$residuals[!wzero]
  Xd <- diag(c(w[!wzero])) %*% ellipticalfit$Xmodel[!wzero, 
                                                    ]
  dev <- 2 * (family$g0(resid, df = family$df, r = family$r, 
                        s = family$s, alpha = family$alpha, mp = family$mp, 
                        epsi = family$epsi, sigmap = family$sigmap, k = family$k) - 
                family$g0(0, df = family$df, r = family$r, s = family$s, 
                          alpha = family$alpha, mp = family$mp, epsi = family$epsi, 
                          sigmap = family$sigmap, k = family$k))
  p <- ellipticalfit$rank
  H <- Xd %*% solve(t(Xd) %*% Wi %*% Xd) %*% t(Xd) %*% Wi
  Hi <- Xd %*% solve(t(Xd) %*% Wi %*% Xd) %*% t(Xd) %*% Wi
  h <- diag(H)/(scalevariance * scale)
  rs <- resid/sqrt(scalevariance * (1 - h))
  ro <- ellipticalfit$y[!wzero] - ellipticalfit$fitted[!wzero]
  n <- length(resid)
  Ones <- rep(1,n)
  u <- ro^2/dispersion
  
  rdesvio <- rep(0,n)
  for(i in 1:n){
    logg0 <- family$g0(0, df = family$df, s = family$s, 
                       r = family$r, alpha = family$alpha, mp = family$mp, 
                       epsi = family$epsi, sigmap = family$sigmap, k = family$k)
    loggu <- family$g0(u[i], df = family$df, s = family$s, 
                       r = family$r, alpha = family$alpha, mp = family$mp, 
                       epsi = family$epsi, sigmap = family$sigmap, k = family$k)
    rdesvio[i] <- sqrt(Wi[i,i])*(sign(resid[i]))*(2*logg0 - 2*loggu)^(.5)
  }
  
  ct <- ellipticalfit$Wgder[!wzero]
  op <- ellipticalfit$Wg[!wzero] * ro
  a <- ellipticalfit$v[!wzero] - 4 * ct * u
  b <- op + (u * ct * ro)
  db <- diag(b)
  b <- matrix(b, n, 1)
  u <- matrix(u, n, 1)
  da <- diag(a)
  dai <- diag(1/a)
  ro <- matrix(ro, n, 1)
  som <- matrix(0, p, p)
  M <- as.matrix(som + t(Xd) %*% da %*% Wi %*% Xd)
  lbb <- (-1/dispersion) * M
  lby <- (1/dispersion) * t(Xd) %*% Wi %*% da
  lphiy <- (-2/(dispersion^2)) * t(b) %*% Wi
  lbphi <- (2/dispersion^2) * t(Xd) %*% Wi %*% b
  lphib <- t(lbphi)
  lphi <- (1/(dispersion^2)) * ((t(Ones)%*%Wi%*%Ones)/2 + t(u) %*% diag(ct) %*% 
                                  Wi %*% u - (1/dispersion) * t(ro) %*% diag(ellipticalfit$v[!wzero]) %*% 
                                  Wi %*% ro)
  lbb1 <- solve(lbb)
  lc1 <- lbb1
  E <- as.vector(lphi - lphib %*% lbb1 %*% lbphi)
  Fi <- -lbb1 %*% lbphi
  GLbeta <- Xd %*% (-lbb1) %*% lby
  R <- Xd %*% (Fi %*% t(Fi) %*% lby + Fi %*% lphiy)
  GLphi <- (-1/E) * R
  G <- GLphir <- matrix(0, n, n)
  
  Om <- rbind(cbind(lbb, lbphi), cbind(lphib, lphi))
  Fr <- -lc1 %*% lbphi
  lc1 <- lc1 + (1/E) * Fr %*% t(Fr)
  lc2 <- (1/E) * Fr
  lc3 <- t(lc2)
  lc4 <- matrix(1/E, 1, 1)
  lc <- cbind(rbind(lc1, lc3), rbind(lc2, lc4))
  GLbeta <- diag(GLbeta)
  GLphi <- diag(GLphi)
  GLphir <- diag(GLphir)
  G <- diag(G)
  GL <- GLbeta + G + GLphi + GLphir
  Bi <- (a/dispersion) * GL
  
  DGbeta <- DGphi <- matrix(0,n,1)
  vu2 <- ellipticalfit$v[!wzero]^2
  for(i in 1:n){
    Delta <- diag(0,n) ; Delta[i,i] = 1
    Delta <- diag(1,n) - Delta; 
    DGbeta[i] <- (vu2[i]/scale*(1 - H[i,i])^2)* u[i]*Wi[i,i]*H[i,i]
    DGphi[i] <- (sum(diag(Wi))*(sqrt(vu2[i])*u[i]-1)^2*Wi[i,i]^2)/((sum(diag(Wi%*%Delta)))^2*scaledispersion)
  }
  
  deltab <- matrix(0, n, p)
  deltad <- matrix(0, n, 1)
  deltab <- (1/dispersion) * diag((ellipticalfit$y[!wzero] - 
                                     ellipticalfit$fitted[!wzero]) * ellipticalfit$v[!wzero]) %*% 
    Wi %*% Xd
  deltad <- Wi %*% matrix(-(0.5/dispersion) * (1 - ellipticalfit$v[!wzero] * 
                                                 u), n, 1)
  delta <- t(cbind(deltab, deltad))
  b11 <- cbind(matrix(0, p, p), matrix(0, p, 1))
  b12 <- cbind(matrix(0, 1, p), 1/E)
  b1 <- rbind(b11, b12)
  b211 <- cbind(lbb1, matrix(0, p, 1))
  b212 <- cbind(matrix(0, 1, p), matrix(0, 1, 1))
  b2 <- rbind(b211, b212)
  Cic <- -t(delta) %*% (lc) %*% delta
  Cic <- 2 * diag(Cic)
  A <- as.matrix(t(delta) %*% (lc) %*% delta)
  decA <- eigen(A)
  Lmax <- decA$val[1]
  dmax <- decA$vec[, 1]
  dmax <- dmax/sqrt(Lmax)
  dmaxc <- abs(dmax)
  
  deltab <- (-2/dispersion) * t(Xd) %*% Wi %*% db
  deltad <- (-1/(dispersion^2)) * t(ro) %*% Wi %*% db
  delta <- rbind(deltab, deltad)
  b11 <- cbind(matrix(0, p, p), matrix(0, p, 1))
  b12 <- cbind(matrix(0, 1, p), 1/lphi)
  b1 <- rbind(b11, b12)
  b211 <- cbind(lbb1, matrix(0, p, 1))
  b212 <- cbind(matrix(0, 1, p), 0)
  b2 <- rbind(b211, b212)
  Ci <- -t(delta) %*% (lc) %*% delta
  Cih <- 2 * diag(Ci)
  A <- as.matrix(t(delta) %*% (lc) %*% delta)
  decA <- eigen(A)
  Lmax <- decA$val[1]
  dmax <- decA$vec[, 1]
  dmax <- dmax/sqrt(Lmax)
  dmax <- abs(dmax)
  
  deltab <- (1/dispersion) * t(Xd) %*% Wi %*% da
  deltad <- (-2/(dispersion^2)) * t(b) %*% Wi
  delta <- rbind(deltab, deltad)
  b11 <- cbind(matrix(0, p, p), matrix(0, p, 1))
  b12 <- cbind(matrix(0, 1, p), 1/lphi)
  b1 <- rbind(b11, b12)
  b211 <- cbind(lbb1, matrix(0, p, 1))
  b212 <- cbind(matrix(0, 1, p), 0)
  b2 <- rbind(b211, b212)
  Ci <- -t(delta) %*% (lc) %*% delta
  Ci <- 2 * diag(Ci)
  ds <- diag(sqrt(dispersion), n)
  deltai <- (1/dispersion) * t(Xd) %*% da %*% Wi %*% ds
  Lmax <- NULL
  A <- matrix(0, n, n)
  for (i in 1:n) {
    A[, i] <- as.matrix(t(deltai) %*% solve(t(Xd) %*% da %*% 
                                              Wi %*% Xd) %*% matrix(Xd[i, ], p, 1))
  }	
  Lmaxr <- abs(diag(A))
  
  Cmaxc <- Lmaxc <- matrix(0, p, n)
  for (j in 1:p) {
    Ff <- matrix(0, p, n)
    Ff[j, ] <- rep(1, n)
    De <- diag(as.vector(ro), n)
    Dv <- diag(ellipticalfit$v[!wzero])
    st <- sqrt(var(Xd[, j]))
    for (i in 1:n) {
      A[, i] <- st * (1/dispersion) * t(Ff %*% De %*% 
                                          Wi %*% Dv - ellipticalfit$coef[j] * t(Xd) %*% Wi %*% da) %*% 
        solve(t(Xd) %*% da %*% Wi %*% Xd) %*% matrix(Xd[i, 
                                                        ], p, 1)
      Cmaxc[j, i] <- 2 * abs(t(A[, i]) %*% A[, i])
    }
    Lmaxc[j, ] <- abs(diag(A))
  }
  
  list(resid = resid, rs = rs, rd = rdesvio, dispersion = dispersion, h = diag(H), Hat = Hi, 
       GL = GL, GLbeta = GLbeta, GLphi = GLphi, DGbeta = DGbeta, DGphi = DGphi, Cih = Cih, 
       Cic = Cic, Lmaxr = Lmaxr, Lmaxc = Lmaxc, Cmaxc = Cmaxc)
}
