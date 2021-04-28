#' @title Multiscale Geographically Weighted Elliptical Regression
#' @import spgwr
#' @import graphics
#' @import methods
#' @import sp
#' @import spData
#' @import maptools
#' @importFrom GWmodel gw.dist gw.weight Generate.formula
#' @name gwer.multiscale
#' @description The function fit geographically weighted elliptical regression model to explore the non-stationarity relationshps across differente spatial scales.
#' @param formula regression model formula as in \code{glm}.
#' @param kernel function chosen as follows:
#' gaussian: wgt = exp(-.5*(vdist/bw)^2);
#' exponential: wgt = exp(-vdist/bw);
#' bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
#' tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
#' boxcar: wgt=1 if dist < bw, wgt=0 otherwise
#' @param adaptive defines the type of bandwidth used. either NULL (default) or a proportion between 0 and 1 of observations to include in weighting scheme (k-nearest neighbours).
#' @param criterion criterion for determining the convergence of the back-fitting procedure, could be "CVR" or "dCVR", which corespond to the changing value of RSS (CVR) and the differential version (dCVR), respectively; and "dCVR" is used as default.
#' @param hatmatrix if TRUE, return the hatmatrix as a component of the result.
#' @param spdisp if TRUE dispersion parameter varies geographically.
#' @param family a description of the error distribution to be used in the model (see \code{\link{family.elliptical}} for details of family functions).
#' @param data model data frame, or may be a SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}.
#' @param dispersion an optional fixed value for dispersion parameter.
#' @param weights an optional numeric vector of weights to be used in the fitting process.
#' @param subset an optional numeric vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain NAs (see \code{glm}).
#' @param longlat TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers. If x is a SpatialPoints object, the value is taken from the object itself.
#' @param control a list of parameters for controlling the fitting process. For \code{elliptical} this is passed by \code{glm.control}.
#' @param max.iterations maximum number of iterations in the back-fitting procedure.
#' @param threshold threshold value to terminate the back-fitting iteration.
#' @param dMats a list of distance matrices used for estimating each specific parameter
#' @param p.vals a collection of positive numbers used as the power of the Minkowski distance
#' @param theta.vals a collection of values used as angles in radians to rotate the coordinate system
#' @param bws0 a vector of initializing bandwidths for the back-fitting procedure, of which the length should equal to the number of paramters if specified
#' @param bw.seled a vector of boolean variables to determine whether the corresponding bandwidth should be re-selected or not: if TRUE, the corresponding bandwiths for the specific parameters are supposed to be given in bws0; otherwise, the bandwidths for the specific parameters will be selected within the back-fitting iterations.
#' @param approach specified by CV for cross-validation approach or by AIC corrected (AICc) approach
#' @param bws.thresholds threshold values to define whether the bandwidth for a specific parameter has converged or not
#' @param bws.reOpts the number times of continually optimizing each parameter-specific bandwidth even though it meets the criterion of convergence, for avoiding sub-optimal choice due to illusion of convergence;
#' @param verbose if TRUE (default) reports the progress of search for bandwidth.
#' @param predictor.centered a logical vector of length equalling to the number of predictors, and note intercept is not included; if the element is TRUE, the corresponding predictor will be centered.
#' @param nlower the minmum number of nearest neighbours if an adaptive kernel is used
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
#' \item{gweight}{the \code{gweights} object used.}
#' \item{spdisp}{the \code{spdisp} object used.}
#' \item{this.call}{the function call used.}
#' \item{fp.given}{the \code{fp.given} object used.}
#' \item{longlat}{the \code{longlat} object used.}
#' @references Brunsdon, C., Fotheringham, A. S. and Charlton, M. E. (1996). 
#' Geographically weighted regression: a method for exploring spatial nonstationarity.
#' Geographical analysis, 28(4), 281-298. \doi{10.1111/j.1538-4632.1996.tb00936.x}
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
#' msgwr.fit.t <- gwer.multiscale(fit.formula, family = Student(3), data = gSRDF, 
#'                                bws0 = rep(gwer.bw.t, 5), hatmatrix = TRUE, 
#'                                adaptive = TRUE)
#' }
#' @rdname gwer.multiscale
#' @export
    
gwer.multiscale <- function (formula, data, kernel = "bisquare", approach = "CV", adaptive = FALSE, criterion = "dCVR", 
                            family = Normal, threshold = 0.00001, dMats, p.vals, theta.vals, longlat = NULL, bws0, 
                            bw.seled = rep(F, length(bws0)), bws.thresholds = rep(0.1, length(dMats)), bws.reOpts = 5, 
                            spdisp = 'local', verbose = F, weights, dispersion = NULL, na.action = "na.fail", hatmatrix = T, 
                            control = glm.control(epsilon = 1e-04, maxit = 100, trace = F), model = FALSE, x = FALSE, 
                            y = TRUE, contrasts = NULL, parplot = FALSE, max.iterations = 2000, subset, offset, 
                            predictor.centered = rep(T, length(bws0) - 1), nlower = 10, ...)
{
  timings <- list()
  timings[["start"]] <- Sys.time()  
  call <- match.call()
  this.call <- match.call()
  dist <- as.character(call$family)[1]
  user.def <- F
  
  if (!charmatch(dist, c("Normal", "Cauchy", "Student", "Gstudent", "LogisI", "LogisII", "Glogis", "Cnormal", "Powerexp"), nomatch = F)) 
    dist <- match.arg(dist, c("Normal", "Cauchy", "Student", "Gstudent", "LogisI", "LogisII", "Glogis", "Cnormal", "Powerexp"))
  else user.def <- T
  
  p4s <- as.character(NA)
  if (is(data, "Spatial")) {
    p4s <- proj4string(data)
    dp.locat <- coords <- coordinates(data)
    regression.points <- data
    data <- as(data, "data.frame")
  }  
  else {
    stop("Given regression data must be Spatial*DataFrame")
  }
  
  if (is.null(longlat) || !is.logical(longlat)) 
    longlat <- FALSE
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
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$family <- mf$control <- mf$model <- mf$dispersion <- mf$x <- mf$y <- mf$kernel <- mf$adaptive <- mf$longlat <-mf$contrasts <- mf$offset <- mf$hatmatrix <- mf$spdisp <- mf$parplot <- mf$... <- NULL
  mf$drop.unused.levels <- TRUE  
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  dp.n <- nrow(dp.locat)  

  if (!missing(family) && !charmatch(dist, c("Normal", "Cauchy", "Student", "Gstudent", "LogisI", "LogisII", "Glogis", "Cnormal", "Powerexp"), F)) 
    cat(paste("\n work with user-defined family:", call$family, "\n"))

  if (!missing(dispersion) && is.numeric(dispersion) && !(dispersion > 0)) 
    stop("\n no negative values for dispersion parameter")


  Y <- model.extract(mf, "response")
  if (!is.numeric(Y)) 
    stop("\n response must be numeric")
  X <- model.matrix(mt, mf, contrasts)
  if (!is.numeric(X)) 
    stop("\n model matrix must be numeric")
  offset <- model.extract(mf, offset)
  nobs <- nrow(X)
  var.n <- ncol(X)

  if (length(offset) == 1 && offset == 0) 
    offset <- rep(0, nobs)
  idx1 <- match("(Intercept)", colnames(X))
  
  w <- model.extract(mf, weights)
  wzero <- rep(F, nrow(mf))

  if (!length(w)) 
    w <- rep(1, nrow(mf))
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
  
#  method <- "gwer.ms.fit" 
#  elliptical.fitter <- get(method)
  offset4fit <- offset

#  fit <- elliptical.fitter(X = X, Y = Y, offset = offset4fit, family = family, 
#                           dispersion = dispersion, maxit = control$maxit, 
#                           epsilon = control$epsilon, trace = control$trace, ...)
  fit <- elliptical(formula =  formula, family = family, data = data, dispersion = dispersion,
                    control = control, offset = offset4fit)
  if (spdisp == 'global')
    dispersion <- fit$dispersion
  
  dispersion <- rep(dispersion, dp.n)
  fit$x <- X
  fit$y <- Y
  
  if (!is.na(idx1)) {
    if (length(predictor.centered) != var.n - 1) {
      predictor.centered <- rep(T, length.out = var.n - 1)
      warning("All the predictors will be centered, please check the parameter predictor.centered")
    }
  }
  else {
    if (length(predictor.centered) != var.n) {
      predictor.centered <- rep(T, length.out = var.n)
      warning("All the predictors will be centered, please check the parameter predictor.centered")
    }
  }
  
  n.cent <- length(which(predictor.centered))
  if(!is.na(idx1)) {
    colnames(X)[idx1] <- "Intercept"
    X1 <- X[,-idx1]
  }
  else {
    X1 <- X
  }
  
  if(n.cent > 1)
    predictors.centered.means <- colMeans(X1[,predictor.centered])
  else if (n.cent == 1)
    predictors.centered.means <- mean(X1[,predictor.centered])
  else 
    predictors.centered.means <- NA
  if(n.cent >= 1) {
    X1[,predictor.centered] <- scale(X1[,predictor.centered], scale=FALSE)
  }
  if(!is.na(idx1))
  {
    X1 <- cbind(1, X1)
  }
  colnames(X1) <-  colnames(X)
  
  allvars <- all.vars(formula)
  DeVar <- allvars[1]
  InDevars <- colnames(X)

  if(missing(dMats)) {
    dMats <- list()
    if(missing(p.vals)) {
      p.vals <- 2
      dMat <- gw.dist(dp.locat, longlat = longlat)
      for(i in 1:var.n)
        dMats[[i]] <- dMat
    }
    else {
      p.vals <- rep_len(p.vals, var.n)
      if(missing(theta.vals))
        theta.vals <- rep_len(0, var.n)
      else
        theta.vals <- rep_len(theta.vals, var.n)
      for(i in 1:var.n){
        dMats[[i]] <- gw.dist(dp.locat, p = p.vals[i], theta = theta.vals[i], longlat = longlat)
      }
    }
  }
  else if(is.list(dMats)) {
    if(length(dMats) != var.n)
      stop("Please specify a distance matrix for each independent variable!")
    else {
      for(i in 1:var.n){
        dim.dMat <- dim(dMats[[i]])
        if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n)
          stop("Dimensions of dMat are not correct")
      }
    }
  }
  else {
    stop("dMats are not correct")
  }  

  
  if(missing(bws0)) {
    bws0 <- numeric(var.n)
    cat("------ Calculate the initial bandwidths for each independent variable ------\n")
    for(i in 1:var.n) {
      if(InDevars[i] == "Intercept")
        fml <- Generate.formula(DeVar, c("1"))
      else {
        fml <- Generate.formula(DeVar, c(InDevars[i]))
      }
      cat("Now select an optimum bandwidth for the model: ", fml, "\n")
      part1 <- paste("bws0[i]<-bw.gwr(", fml, sep="")
      part2 <- "data=regression.points,kernel=kernel,approach=approach,dMat=dMats[[i]])"
      expression <- paste(part1, part2, sep=",")
      print(expression)
      eval(parse(text = expression))
#      bws0[i] <- bw.gwer(X1, Y, dp.locat, approach = approach, kernel = kernel, adaptive = adaptive, 
#                        dMat, verbose = verbose, nlower = nlower)
    }
    cat("------            The end for the initial selections              ------\n")
  }
  else {
    bws0 <- rep(bws0, length.out = var.n)
  }
  

  if(length(bw.seled) != var.n)
    bw.seled <- rep(F, length.out = var.n)
  cat("------   Calculate the initial beta0 from the above bandwidths    ------\n")
  dMat <- dMats[[1]]
  bw.int0 <- bw.gwer1(X1, Y, dp.locat = dp.locat, approach = approach, kernel = kernel, adaptive = adaptive, 
                      dispersion = dispersion, dMat = dMat, family = family, verbose = verbose, nlower = nlower)

  if(hatmatrix) {
    Shat <- matrix(nrow = dp.n, ncol = dp.n)
    S.arrays <- array(dim = c(var.n, dp.n, dp.n))
    H <- matrix(nrow = dp.n, ncol = dp.n)
    H.arrays <- array(dim = c(var.n, dp.n, dp.n))
    C <- array(dim=c(dp.n, var.n, dp.n))
  }
  ms.gweights <- array(dim = c(var.n, dp.n, dp.n))

  fit.gwer <- gwer.mfit(X1, Y, family = family, loc = dp.locat, adaptive = adaptive, hatmatrix = hatmatrix, bw = bw.int0, kernel = kernel, 
                   spdisp = spdisp, dMat = dMat, offset = offset4fit, dispersion = dispersion, control = control)
  betas <- fit.gwer[[1]]
  #std.error <- fit.gwer[[8]]
  #p.values <- fit.gwer[[5]]

  if(spdisp == 'local' || spdisp == 'global')
    dispersion <- fit.gwer[[4]]
  
  if(hatmatrix) {
    Shat <- fit.gwer[[2]]
    C <- fit.gwer[[3]]
    H <- fit.gwer[[6]]
    Haux <- fit.gwer[[7]]
    idm <- diag(var.n)
    for(i in 1:var.n)
      for(j in 1:dp.n)
      {
        S.arrays[i, j, ] <- X1[j, i]*(idm[i, ]%*%C[j, , ])
        H.arrays[i, j, ] <- X1[j, i]*(idm[i, ]%*%Haux[j, , ])
      }
  }

  cat("------            The end for calculating the initial beta0              ------\n") 
  cat("------ Select the optimum bandwidths for each independent variable via ", approach, " aproach ------\n")
  iteration <-0 
  bws.vars <- bws0
  bws.change.NO <- numeric(var.n)
  criterion.val <- 10000000  
  resid.i <- Y - (X1*betas)%*%matrix(1, ncol(X1), 1)#####gw.fitted(X1, betas)
  RSS0 <- sum(resid.i^2)
  RSS1 <- 0
  RSS.vals <- c(RSS0, RSS1, criterion.val)
  cat("*****  The back-fitting process for model calibration with bandwiths selected *****\n")
  
#  gweights <- matrix(0,n,n) ; local.fitted <- matrix(0,n,n) ; sum.w <- numeric(n)# ; resid.ord <- matrix(0,n,1)
#  
#  if(spdisp){ 
#    dispersion.gwer <- numeric(n) 
#    gwr.b <- matrix(nrow = nobs, ncol = ncol(X)+1)
#    stderror <- matrix(0,n,p+1) ; zvalue <- matrix(0,n,p) ; pvalue <- matrix(0,n,p)
#    colnames(gwr.b) <- c(colnames(X), "dispersion")
#    colnames(stderror) <- c(paste(colnames(X),rep(".se",p), sep = ""), "dispersion.se")
#  } else {
#    gwr.b <- matrix(nrow = nobs, ncol = ncol(X))
#    stderror <- matrix(0,n,p) ; zvalue <- matrix(0,n,p) ; pvalue <- matrix(0,n,p)
#    colnames(gwr.b) <- c(colnames(X))    
#    colnames(stderror) <- paste(colnames(X),rep(".se",p), sep = "")
#  }
#  colnames(zvalue) <- colnames(pvalue) <- c(colnames(X)) 

  # Iterative process #
  while((iteration < max.iterations) && criterion.val > threshold) {
    cat("    Iteration ", iteration + 1, ":\n")   
    for(i in 1:var.n) {
      dMat <- dMats[[i]]  
      f.i <- betas[, i]*X1[, i]
      X.i <- matrix(X1[, i], ncol = 1)
      colnames(X.i) <- colnames(X1)[i]
      Y.i <- resid.i + f.i
      if(bw.seled[i])
        bw.i <- bws0[i]
      else {
        cat("Now select an optimum bandwidth for the variable: ", InDevars[i], "\n")

        bw.i <- bw.gwer1(X.i, Y.i, dp.locat = dp.locat, approach = approach, kernel = kernel, adaptive = adaptive, 
                         dispersion = dispersion, dMat = dMat, family = family, verbose = verbose, nlower = nlower)
        cat("The newly selected bandwidth for variable ", InDevars[i])
        cat(" is: ", bw.i, "\n")
        cat("The bandwidth used in the last iteration is: ", bws0[i])
        cat(" and the difference between these two bandwidths is: ", abs(bw.i - bws0[i]), "\n")
        if (abs(bw.i - bws0[i]) > bws.thresholds[i]) {
          cat("The bandwidth for variable ", InDevars[i])
          cat(" will be continually selected in the next iteration.\n")
          bws.change.NO[i] <- 0
        }
        else {
          bws.change.NO[i] <- bws.change.NO[i] + 1
          if(bws.change.NO[i] < bws.reOpts) {
            cat("The bandwidth for variable ", InDevars[i])
            cat(" seems to be converged for ", bws.change.NO[i], " times.")
            cat("It will be continually optimized in the next ", bws.reOpts-bws.change.NO[i], " times\n")
          }
          else {
            cat("The bandwidth for variable ", InDevars[i])
            cat(" seems to be converged and will be kept the same in the following iterations.\n")
            bw.seled[i] <- T
          }
        }
      }
      bws0[i] <- bw.i       

      #res <- gwr.q2(matrix(X1[,i], ncol=1), y.i, dp.locat, adaptive=adaptive, hatmatrix = hatmatrix,bw=bw.i, kernel=kernel,dMat=dMat)
      fit.i <- gwer.mfit(X.i, Y.i, family = family, loc = dp.locat, adaptive = adaptive, hatmatrix = hatmatrix, 
                        bw = bw.i, kernel = kernel, dMat = dMat, spdisp = spdisp, offset = offset4fit, 
                        dispersion = dispersion, control = control)
      beta.i <- fit.i[[1]]
#      std.error.i <- fit.i[[8]]
#      p.values.i <- fit.i[[5]]
      
#      dxs <- spDistsN1(coords, fit.points[i, ], longlat = longlat)
#      if (any(!is.finite(dxs))) 
#        dxs[which(!is.finite(dxs))] <- .Machine$double.xmax/2
#      w.i <- gweight(dxs^2, bandwidth[i])
#      if (any(w.i < 0 | is.na(w.i))) 
#        stop(paste("Invalid weights for i:", i))
#      fit.i <- elliptical.fitter(X = X, Y = Y, gweights = w.i, offset = offset4fit, 
#                                family = family, dispersion = dispersion, 
#                                maxit = control$maxit, epsilon = control$epsilon, 
#                                trace = control$trace, ...)
#      sum.w[i] <- sum(w.i)
#      v_fitteds[i] <- fit.i$fitted.values[i]
#      if(spdisp){ 
#        dispersion.gwer[i] <- fit.i$dispersion
#        gwr.b[i, ] <- c(coefficients(fit.i), dispersion.gwer[i])
#      } else{
#        gwr.b[i, ] <- coefficients(fit.i)
#      }
#      if (!fp.given && hatmatrix) {   
#        if(spdisp){
#          stderror[i, ] <- c(rowlen[1:p] * sqrt(fit.i$dispersion/fit.i$scale), sqrt(4*fit.i$dispersion^2/(sum(fit.i$gweights)*fit.i$scaledispersion)))
#          zvalue[i, ] <- c(gwr.b[i, 1:p])/stderror[i, 1:p]
#          pvalue[i, ] <- 2 * pnorm(-abs(zvalue[i, ]))
#        } else {
#          stderror[i, ] <- rowlen[1:p] * sqrt(fit.i$dispersion/fit.i$scale)
#          zvalue[i, ] <- c(gwr.b[i, ])/stderror[i, ]
#          pvalue[i, ] <- 2 * pnorm(-abs(zvalue[i, ]))
#        }
#        Hat[i,] <- X[i,] %*% solve(t(X) %*% diag(w.i) %*% X) %*% t(X) %*% diag(w.i)    
#      }
      if(spdisp == 'multi'){
        disp.mult[, i] <- fit.i[[4]]
        dispersion <- apply(disp.mult, 2, mean)
      }
      
      if(hatmatrix)
      {
        Si <- fit.i[[2]]
        S.arrayi <- S.arrays[i,,]
        S.arrays[i,,] <- Si%*%S.arrayi + Si - Si%*%Shat
        Shat <- Shat - S.arrayi + S.arrays[i,,]
        
        Hi <- fit.i[[6]]
        H.arrayi <- H.arrays[i,,]
        H.arrays[i,,] <- Hi%*%H.arrayi + Hi - Hi%*%H
        H <- H - H.arrayi + H.arrays[i,,]
      }
      
      betas[, i] <- beta.i
      #std.error[, i] <- std.error.i
      #p.values[, i] <- p.values.i/(sum(diag(S.arrays[i,,]))/var.n)
      resid.i <- Y - (X1*betas)%*%matrix(1, ncol(X1), 1) 
      #resid.i <- Y - gw.fitted(X1, betas)
      ms.gweights[i,,] <- fit.i[[9]]
    }
    bws.vars <- rbind(bws.vars, bws0)
    RSS1 <- sum((Y - (X1*betas)%*%matrix(1, ncol(X1), 1))^2)
    #RSS1 <- sum((Y - gw.fitted(X1, betas))^2)   
    if(criterion=="CVR") {
      criterion.val <- abs(RSS1-RSS0)
      cat("    Iteration ", iteration, "the change value of RSS (CVR) is: ", criterion.val,"\n")
    }
    else {
      criterion.val <- sqrt(abs(RSS1-RSS0)/RSS1)
      cat("    Iteration ", iteration + 1, "the differential change value of RSS (dCVR) is: ", criterion.val,"\n") 
    }
 
    RSS0 <- RSS1  
    cat("----------End of    Iteration ", iteration + 1, "----------\n") 
    RSS.vals <- rbind(RSS.vals, c(RSS0, RSS1, criterion.val)) 
    iteration <- iteration + 1            
  }  
  
  # Output #
  local.fitted <- betas %*% t(X1)
  yhat <- (X1*betas)%*%matrix(1, ncol(X1), 1)
  residual <- Y - yhat
  dispersion.gwer <- fit.gwer[[4]]
  
  RSS.gw <- RSS1  
  yss.g <- sum((y - mean(y))^2)
  R2.val <-  1-RSS.gw/yss.g
  R2adj <- NA
  AIC <- NA ; AICc <- NA ; BIC <- NA
  if(hatmatrix) {
    scalevariance <- family$g4(res, df = family$df, 
                               r = family$r, s = family$s, alpha = family$alpha, 
                               mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                               k = family$k)
    D.phi <- diag(dispersion.gwer) ; alpha.pv <- rep(1, var.n)
    std.error <- z.value <- p.values <- matrix(nrow=dp.n, ncol=var.n)
    for(i in 1:var.n) {
      C.i <- solve(diag(X1[, i]))%*%S.arrays[i,,]
      var.beta.i <- C.i%*%D.phi%*%t(C.i)*scalevariance
      std.error[, i] <- sqrt(diag(var.beta.i))
      z.value[, i] <- c(betas[, i])/std.error[, i]
      p.values.i <- 2*pnorm(-abs(z.value[, i]))
      p.values[, i] <- p.values.i#/(sum(diag(S.arrays[i,,]))/var.n)
      alpha.pv[i] <- (sum(diag(S.arrays[i,,]))/var.n)
    }
    nu1 <- sum(diag(Shat))
    nu2 <- sum(diag(Shat^2))

    res <- (Y - yhat)/dispersion.gwer
    rpearson <- res/sqrt(scalevariance)
    #H1 <- (1/(scalevariance * scale)) * H
    #varr <- scalevariance * (1 - diag(H1))
    #rstand <- res/sqrt(varr)
    
    if (charmatch(dist, "Normal", F)) {
      dist.y <- pnorm(Y, yhat, dispersion.gwer)
    }
    else if (charmatch(dist, "Cauchy", F)) {
      dist.y <- pcauchy(Y, yhat, dispersion.gwer)
    }
    else if (charmatch(dist, "Student", F)) {
      dist.y <- pt(res, family$df)
    }
    else if (charmatch(dist, "Gstudent", F)) {
      dist.y <- pgstudent(res, family$s, family$r)
    }
    else if (charmatch(dist, "LogisI", F)) {
      stop(paste("not implemented yet"))
      dist.y <- plogisI(Y, yhat, dispersion.gwer)
    }
    else if (charmatch(dist, "LogisII", F)) {
      dist.y <- plogisII(res)
    }
    else if (charmatch(dist, "Glogis", F)) {
      stop(paste("not implement yet"))
      dist.y <- NULL#pglogis(res, family$alpha, family$mp)
    }
    else if (charmatch(dist, "Cnormal", F)) {
      stop(paste("not implemented yet"))
      dist.y <- NULL#pcnormal(res, family$epsi, family$sigmap)
    }
    else if (charmatch(dist, "Powerexp", F)) {
      dist.y <- ppowerexp(res, family$k)
    }
    rquant <- qnorm(dist.y, 0, 1)
    
    logLik <- -0.5 * sum(log(dispersion.gwer)) + sum(family$g0(res, df = family$df, s = family$s, r = family$r,
                                                               alpha = family$alpha, mp = family$mp, epsi = family$epsi,
                                                               sigmap = family$sigmap, k = family$k))
    edf <- dp.n - 2*nu1 + nu2
    AIC <- 2*nu1 - 2*logLik
    AICc <- -2*logLik + (dp.n*(dp.n + nu1))/(dp.n-2-nu1)
    BIC <- log(dp.n)*nu1 - 2*logLik
    R2.adj <- 1-(1-R2.val)*(dp.n-1)/(edf-1)
  }
  
  GW.diagnostic <- list(RSS.gw = RSS.gw, logLik = logLik, AIC = AIC, AICc = AICc, BIC = BIC, Hat = Shat, R2.val = R2.val, R2.adj = R2.adj)

  if(!is.na(idx1) && n.cent>=1) {
    idx.cent <- which(predictor.centered)
    betas.cent <- c()
    for(i in 1:n.cent)
      betas.cent <- cbind(betas.cent, betas[,idx.cent[i]+1]*predictors.centered.means[i])
    beta0 <- betas[,idx1] - apply(betas.cent, 1, sum)
    betas[,idx1] <- beta0
  }
  
  vdgwer.df <- data.frame(betas, yhat, residual)
  colnames(vdgwer.df) <- c(colnames(X), "yhat", "residual")
  griddedObj <- F
  if (is(regression.points, "Spatial"))
  { 
    if (is(regression.points, "SpatialPolygonsDataFrame"))
    {
      polygons<-polygons(regression.points)
      SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=vdgwer.df,match.ID=F)
    }
    else
    {
      griddedObj <- gridded(regression.points)
      SDF <- SpatialPointsDataFrame(coords=dp.locat, data=vdgwer.df, proj4string=CRS(p4s), match.ID=F)
      gridded(SDF) <- griddedObj 
    }
  }
  else
    SDF <- SpatialPointsDataFrame(coords=dp.locat, data=vdgwer.df, proj4string=CRS(p4s), match.ID=F)  
  timings[["stop"]] <- Sys.time()
  GW.arguments <- list(formula = formula, criterion = criterion, bws = bws0, alpha.pv = alpha.pv, kernel = kernel, adaptive = adaptive, hatmatrix = hatmatrix)
  GW.residual <- list(ordinal = residual, response = res, pearson = rpearson, quantile = rquant)
  
  z <- list(SDF = SDF, coef = betas, se = std.error, pv = p.values, dispersion = dispersion.gwer, gweights = ms.gweights, fitted = yhat, flm = local.fitted, lm = fit, 
            bandwidth = bws.vars, GW.residual = GW.residual, GW.arguments = GW.arguments, GW.diagnostic = GW.diagnostic, family = fit$family, timings = timings, 
            hatmatrix = hatmatrix, this.call = this.call)
  class(z) <- "multiscalegwer"
  invisible(z) 
}



#--------------------#
# Estimation Process #
#--------------------#
gwer.mfit <- function(X, Y, family, loc, adaptive = FALSE, hatmatrix, bw=sqrt(var(loc[,1])+var(loc[,2])),
                      kernel, p, theta, longlat, dMat, wt2=rep(1,nrow(loc)), spdisp, offset, dispersion, control)
{
  if (missing(dMat))
    DM.given <- F
  else
    DM.given <- T
  dp.n <- nrow(loc) ; var.n <- ncol(X)
  gweights <- matrix(nrow=dp.n, ncol=dp.n)
  gwr.b <- matrix(nrow=dp.n, ncol=var.n) ; dispersion.gwer <- rep(0, dp.n) 
  pvalue <- zvalue <- stderror <- matrix(nrow=dp.n, ncol=var.n)
  S <- matrix(nrow=dp.n, ncol=dp.n)
  C <- array(dim=c(dp.n, var.n, dp.n))
  H <- matrix(nrow=dp.n, ncol=dp.n)
  Haux <- array(dim=c(dp.n, var.n, dp.n))
  for (i in 1:dp.n)
  {
    if(DM.given)
      dist.vi <- dMat[,i]
    else
      dist.vi <- gw.dist(loc, focus=i, p, theta, longlat)
    w.i <- gw.weight(dist.vi, bw, kernel, adaptive)
    fit.i <- gwer.fit(X = X, Y = Y, gweights = w.i, offset = offset, 
                      family = family, dispersion = dispersion[i], 
                      maxit = control$maxit, epsilon = control$epsilon, 
                      trace = control$trace)
    #gw.resi <- gw_reg(x,y,as.vector(W.i*wt2),hatmatrix=hatmatrix,i)
    #betas[i,] <- gw.resi[[1]]
    #if(spdisp == 'local' || spdisp == 'multiscale'){ 
    #  dispersion.gwer[i] <- fit.i$dispersion
    #  gwr.b[i, ] <- coefficients(fit.i)#c(coefficients(fit.i), dispersion.gwer[i])
    #} else{
    #  dispersion.gwer[i] <- fit.i$dispersion
    #  gwr.b[i, ] <- coefficients(fit.i)
    #}
    dispersion.gwer[i] <- fit.i$dispersion
    gwr.b[i, ] <- coefficients(fit.i)
    gweights[i,] <- fit.i$gweights
    
    R <- fit.i$R[(1:var.n), (1:var.n)]
    covun <- solve(qr(R))
    rowlen <- sqrt(diag(covun))
    
    if(hatmatrix){
      stderror[i, ] <- rowlen[1:var.n] * sqrt(fit.i$dispersion/fit.i$scale)
      zvalue[i, ] <- c(gwr.b[i, ])/stderror[i, ]
      pvalue[i, ] <- 2 * pnorm(-abs(zvalue[i, ]))
      S[i,] <- fit.i$Hat[i, ]
      C[i,,] <- fit.i$Ci
      H[i,] <- fit.i$H[i, ]
      Haux[i,,] <- fit.i$Haux
      #S[i,]<-gw.resi[[2]]
      #C[i,,]<-gw.resi[[3]]
    }
    #fitteds[i] <- fit.i$fitted.values[i]
  }
  colnames(gwr.b) <- colnames(X)

  fit <- list(gwr.b, S, C, dispersion.gwer, pvalue, H, Haux, stderror, gweights)
  fit
}



#---------------------#
# Bandwidth Selection #
#---------------------#
bw.gwer1 <- function(x, y, dp.locat, approach="AIC", kernel="gaussian", adaptive = FALSE, dispersion, dMat, family, verbose=F, nlower = 10)
{
  dp.n <-  nrow(dp.locat)
  if(adaptive)
  {
    upper <- dp.n
    lower <- nlower
  }
  else
  {
    upper<-range(dMat)[2]
    lower<-upper/5000
  }
  
  #Select the bandwidth by golden selection
  bw<-NA
  if (approach=="cv"||approach=="CV")
    bw <- gold2(gwer.cv, lower, upper, adapt.bw=adaptive, x, y, kernel, adaptive, dispersion, dp.locat, dMat=dMat, family = family, verbose=verbose)
  else if (approach=="aic"||approach=="AIC"||approach=="AICc")
    bw <- gold2(gwer.aic, lower, upper, adapt.bw=adaptive, x, y, kernel, adaptive, dispersion, dp.locat, dMat=dMat, family = family, verbose=verbose)
  else if (approach=="mi"||approach=="MI")
    bw <- gold2(gwer.mi, lower, upper, adapt.bw=adaptive, x, y, kernel, adaptive, dispersion, dp.locat, dMat=dMat, family = family, verbose=verbose)
  bw
}



 
  
