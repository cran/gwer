#' @title Monte Carlo (randomisation) Test for Significance of GWER Parameter Variability
#' @importFrom GWmodel gw.dist gw.weight
#' @name gwer.montecarlo
#' @description This function implements a Monte Carlo (randomisation) test to test for significant (spatial) variability of a geographically weighted elliptical regression model's parameters or coefficients.
#' @param formula regression model formula of a formula \code{object}.
#' @param family a description of the error distribution to be used in the model (see \code{\link{family.elliptical}} for details of elliptical distribution).
#' @param data an optional data frame, list or environment containing the variables in the model.
#' @param nsims the number of randomisations.
#' @param dispersion an optional fixed value for dispersion parameter.
#' @param kernel function chosen as follows:
#' gaussian: wgt = exp(-.5*(vdist/bw)^2);
#' exponential: wgt = exp(-vdist/bw);
#' bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
#' tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
#' boxcar: wgt=1 if dist < bw, wgt=0 otherwise.
#' @param adaptive if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance).
#' @param bw value of the selected bandwidth used in the weighting function (see \code{\link{bw.gwer}} for bandwidth optimization).
#' @param p	 the power of the Minkowski distance, default is 2 (Euclidean distance).
#' @param theta an angle in radians to rotate the coordinate system, default is 0
#' @param longlat if TRUE, great circle distances will be calculated.
#' @param dMat a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}.
#' @param control a list of parameters for controlling the fitting process. This is passed by \code{\link{glm.control}}.
#' @return A vector containing p-values for all parameters spatial variability tests
#' @references Brunsdon C, Fotheringham AS, Charlton ME (1998) Geographically weighted regression - modelling spatial non-stationarity. Journal of the Royal
#' Statistical Society, Series D-The Statistician 47(3):431-443
#' @seealso \code{\link{bw.gwer}}, \code{\link{elliptical}}, \code{\link{family.elliptical}}
#' @keywords Geographically weighted regression
#' @keywords Parameter Variability Tests
#' @keywords Elliptical model
#' @examples
#' \donttest{
#' data(georgia, package = "spgwr")
#' fit.formula <- PctBach ~ TotPop90 + PctRural + PctFB + PctPov
#' gwer.bw.t <- bw.gwer(fit.formula, data = gSRDF, family = Student(4), adapt = TRUE)
#' gwer.fit.t <- gwer(fit.formula, data = gSRDF, family = Student(4), bandwidth = gwer.bw.t, 
#'                    adapt = TRUE, parplot = FALSE, hatmatrix = TRUE, spdisp = TRUE, 
#'                    method = "gwer.fit")
#' gwer.montecarlo(fit.formula, data = gSRDF, family = Student(3), bw = gwer.bw.t, adaptive = TRUE) 
#' }
#' @rdname gwer.montecarlo
#' @export
#' 
gwer.montecarlo <- function (formula, family = Normal, data = list(), nsims = 99, kernel = "bisquare", 
                            adaptive = F, bw, p = 2, theta = 0, dispersion = NULL, longlat = F, dMat, 
                            control = glm.control(epsilon = 1e-04, maxit = 100, trace = F)) 
{
  this.call <- match.call()
  if (!is.null(data)) {
    if (is(data, "Spatial")) {
      dp.locat <- coordinates(data)
      data <- as(data, "data.frame")
      dp.n <- nrow(dp.locat)
    }
    else {
      if (!is(data, "data.frame")) 
        stop("Given regression data must be data.frame or Spatial*DataFrame")
    }
  }
  else stop("No regression data frame is avaiable!")
  if (missing(dMat)) {
    dMat <- gw.dist(dp.locat = dp.locat, p = p, theta = theta, 
                    longlat = longlat)
  }
  else {
    dim.dMat <- dim(dMat)
    if (!is.numeric(dMat) || !is(dMat, "matrix")) 
      stop("Distance matrix(dMat) has to be specified correctly")
    if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n) 
      stop("Dimensions of dMat are not correct")
  }
  if (missing(data)) 
    data <- environment(formula)  
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  dp.n <- length(model.extract(mf, "response"))
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  var.n <- ncol(x)
  if (missing(bw)) 
    stop("Bandwidth must be given for non-adaptive weights")
  if (adaptive) {
    stopifnot(is.numeric(bw))
    stopifnot((bw > 0))
    stopifnot((bw <= dp.n))
  }
  else {
    stopifnot(is.numeric(bw))
    stopifnot((bw > min(dMat)))
  }
  bandwidth <- bw
  betas <- matrix(nrow = dp.n, ncol = var.n)
  Var_betas <- matrix(nrow = nsims + 1, ncol = var.n)
  for (i in 1:dp.n) {
    dist.vi <- dMat[, i]
    w.i <- gw.weight(dist.vi, bw, kernel, adaptive)
    fit.i <- gwer.fit(X = x, Y = y, gweights = w.i, family = family, 
                       dispersion = dispersion, offset = rep(0, dp.n), 
                       maxit = control$maxit, epsilon = control$epsilon, 
                       trace = control$trace)
    betas[i, ] <- coefficients(fit.i)
  }
  for (j in 1:var.n) {
    Var_betas[1, j] <- var(betas[, j])
  }
  sq <- 1:length(y)
  for (k in 1:nsims) {
    mcs <- sample(dp.n)
    dMat[mcs, ] <- dMat[1:dp.n, ]
    dMat[, mcs] <- dMat[, 1:dp.n]
    betas <- matrix(nrow = dp.n, ncol = var.n)
    for (i in 1:dp.n) {
      dist.vi <- dMat[, i]
      w.i <- gw.weight(dist.vi, bw, kernel, adaptive)
      fit.i <- gwer.fit(X = x, Y = y, gweights = w.i, family = family, 
                         dispersion = dispersion, offset = rep(0, dp.n), 
                         maxit = control$maxit, epsilon = control$epsilon, 
                         trace = control$trace)
      betas[i, ] <- coefficients(fit.i)
    }
    for (j in 1:var.n) {
      Var_betas[k + 1, j] <- var(betas[, j])
    }
  }
  p.values <- numeric(var.n)
  for (j in 1:var.n) {
    var_betaj <- Var_betas[, j]
    indx <- sort(var_betaj, index.return = T)$ix
    p.values[j] = 1 - (which(indx == 1))/(nsims + 1)
  }
  pmat <- matrix(p.values, ncol = 1)
  colnames(pmat) <- c("p-value")
  rownames(pmat) <- colnames(x)
  cat("\nTests based on the Monte Carlo significance test\n\n")
  printCoefmat(pmat)
  invisible(pmat)
} 