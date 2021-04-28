#' @title Diagnostic for Geographically Weighted Elliptical Regression Models
#' @description This function obtains the values of different residuals types and calculates the diagnostic measures for the fitted geographically weighted elliptical regression model.
#' @param object an object with the result of the fitted multiscale geographically weighted elliptical regression model.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return Returns a list of diagnostic arrays:
#' \item{ro}{ordinal residuals.}
#' \item{rr}{response residuals.}
#' \item{rp}{pearson residuals.}
#' \item{rs}{studentized residuals.}
#' \item{rd}{deviance residuals.}
#' \item{dispersion}{coefficient of dispersion parameter.}
#' \item{Hat}{the hat matrix.}
#' \item{h}{main diagonal of the hat matrix.}
#' \item{GL}{generalized leverage.}
#' \item{GLbeta}{generalized leverage of location parameters estimation.}
#' \item{GLphi}{generalized leverage of dispersion parameters estimation.}  
#' \item{DGbeta}{cook distance of location parameters estimation.}
#' \item{DGphi}{cook distance of dispersion parameters estimation.}  
#' \item{Cic}{normal curvature for case-weight perturbation.}
#' \item{Cih}{normal curvature for scale perturbation.}
#' \item{Lmaxr}{local influence on response (additive perturbation in responce).}
#' \item{Lmaxc}{local influence on coefficients (additive perturbation in predictors).}
#' @references Brunsdon, C., Fotheringham, A. S. and Charlton, M. E. (1996). 
#' Geographically weighted regression: a method for exploring spatial nonstationarity.
#' Geographical analysis, 28(4), 281-298. \doi{10.1111/j.1538-4632.1996.tb00936.x}
#' @references Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
#' symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
#' \doi{10.1016/j.spl.2005.04.033}
#' @seealso \code{\link{elliptical}}
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
#' gwer.multiscale.diag(msgwr.fit.t)
#' }
#' @export

gwer.multiscale.diag <- function (object,...) 
{
  if (!object$hatmatrix) 
    stop("Diagnostic measures not applicable - regression points different from observed points")
  
  Wik <- object$gweights ; l.fitted = object$flm
  pk <- object$lm$rank ; p <- 1
  family <- object$family
  dist <- family[[1]]
  user.def <- object$lm$user.def
  f.name <- family[[1]]
  dispersion <- as.numeric(object$dispersion)
  w <- if (is.null(object$lm$weights)) 
    rep(1, length(object$lm$residuals))
  else object$lm$weights
  wzero <- (w == 0)
  Xdk <- diag(c(w[!wzero])) %*% object$lm$Xmodel[!wzero, 
                                                       ]
  n <- length(object$lm$residuals)
  if(length(dispersion)==1)
    dispersion <- rep(dispersion, n)
  ro <- rr <- rp <- rq <- rep(0,n) ; H <- matrix(0,n,n)#rp <- rs <- rd <- rep(0,n)
  #H <- G <- GL <- GLbeta <- GLphi <- GLphir <- Bi <- Lmaxr <- Cic <- Cih <- DGbeta <- DGphi <- DGphir <- matrix(0,n,n)
  #Lmaxr <- Cic <- Cih <- array(0, pk, n, n) ; Cmaxc <- Lmaxc <- matrix(0, pk, n)
  Lmaxr <- Cic <- Cih <- Cir <- matrix(0,n,n) ; Cmaxc <- Lmaxc <- matrix(0, pk, n)
  Cic <- list() ; for(j in 1:pk){Cic[[j]] <- matrix(0,n,n)}
  Cih <- list() ; for(j in 1:pk){Cih[[j]] <- matrix(0,n,n)}
  Cir <- list() ; for(j in 1:pk){Cir[[j]] <- matrix(0,n,n)}
  Lmaxr <- list() ; for(j in 1:pk){Lmaxr[[j]] <- matrix(0,n,n)}
  Lmaxc <- list() ; for(j in 1:pk){Lmaxc[[j]] <- matrix(0,n,n)}
  for(k in 1:pk){
    Wi <- Wik[k, , ]
    Xd <- matrix(Xdk[ ,k], n, 1)#as.matrix(Xdk[ ,k])
    for(i in 1:n){
      roi <- object$lm$y[!wzero] - object$flm[i,!wzero]
      resid <- (object$lm$y[!wzero] - object$flm[i,!wzero])/sqrt(dispersion[i])
      Ones <- rep(1,n) ; u <- roi^2/dispersion[i]
      scale <- 4 * family$g2(resid, df = family$df, 
                           r = family$r, s = family$s, alpha = family$alpha, 
                           mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                           k = family$k)    
      scaledispersion <- -1 + 4 * family$g3(resid, 
                                            df = family$df, r = family$r, s = family$s, alpha = family$alpha, 
                                            mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                                            k = family$k)
      scalevariance <- family$g4(resid, df = family$df, 
                                r = family$r, s = family$s, alpha = family$alpha, 
                                mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                                k = family$k)
      Wg <- family$g1(resid, df = family$df, 
                      r = family$r, s = family$s, alpha = family$alpha, 
                     mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                      k = family$k)
      Wgder <- family$g5(resid, df = family$df, 
                        r = family$r, s = family$s, alpha = family$alpha, 
                        mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                        k = family$k)
      V <- -2 * family$g1(resid, df = family$df, 
                          r = family$r, s = family$s, alpha = family$alpha, 
                          mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                          k = family$k)

      ro[i] <- roi[i] ; rr[i] <- resid[i] ; rp[i] <- resid[i]/sqrt(scalevariance)
    
      dev <- 2 * (family$g0(resid, df = family$df, r = family$r, 
                            s = family$s, alpha = family$alpha, mp = family$mp, 
                            epsi = family$epsi, sigmap = family$sigmap, k = family$k) - 
                  family$g0(0, df = family$df, r = family$r, s = family$s, 
                            alpha = family$alpha, mp = family$mp, epsi = family$epsi, 
                            sigmap = family$sigmap, k = family$k))
      wi <- diag(Wi[i,])
      #Hi <- Xd %*% solve(t(Xd) %*% wi %*% Xd) %*% t(Xd) %*% wi
      #H[i,] <- Xd[i,] %*% solve(t(Xd) %*% wi %*% Xd) %*% t(Xd) %*% wi
      #h <- diag(Hi)/(scalevariance * scale)
      #rs[i] <- resid[i]/sqrt(scalevariance * (1 - h[i]))

      #logg0 <- family$g0(0, df = family$df, s = family$s, 
      #                   r = family$r, alpha = family$alpha, mp = family$mp, 
      #                   epsi = family$epsi, sigmap = family$sigmap, k = family$k)
      #loggu <- family$g0(u[i], df = family$df, s = family$s, 
      #                   r = family$r, alpha = family$alpha, mp = family$mp, 
      #                   epsi = family$epsi, sigmap = family$sigmap, k = family$k)
      #rd[i] <- sqrt(Wi[i,i])*(sign(resid[i]))*(2*logg0 - 2*loggu)^(.5)
  

      ct <- Wgder[!wzero]
      op <- Wg[!wzero] * roi
      a <- V[!wzero] - 4 * ct * u
      b <- op + (u * ct * roi)
      db <- diag(b)
      b <- matrix(b, n, 1)
      u <- matrix(u, n, 1)
      da <- diag(a)
      #dai <- diag(1/a)
      roi <- matrix(roi, n, 1)
      som <- matrix(0, p, p)
      #Gi <- GLphiri <- matrix(0, n, n)
      deltab <- matrix(0, n, p)
      deltad <- matrix(0, n, 1)

      M <- as.matrix(som + t(Xd) %*% da %*% wi %*% Xd)
      lbb <- (-1/dispersion[i]) * M
      lby <- (1/dispersion[i]) * t(Xd) %*% wi %*% da
      lphiy <- (-2/(dispersion[i]^2)) * t(b) %*% wi
      lbphi <- (2/dispersion[i]^2) * t(Xd) %*% wi %*% b
      lphib <- t(lbphi)
      lphi <- (1/(dispersion[i]^2)) * ((t(Ones)%*%wi%*%Ones)/2 + t(u) %*% diag(ct) %*% 
                                    wi %*% u - (1/dispersion[i]) * t(roi) %*% diag(V[!wzero]) %*% 
                                    wi %*% roi)
      lbb1 <- solve(lbb)
      lc1 <- lbb1
      E <- as.vector(lphi - lphib %*% lbb1 %*% lbphi)
      #Fi <- -lbb1 %*% lbphi
      #GLbetai <- Xd %*% (-lbb1) %*% lby
      #R <- Xd %*% (Fi %*% t(Fi) %*% lby + Fi %*% lphiy)
      #GLphii <- (-1/E) * R
  
      Om <- rbind(cbind(lbb, lbphi), cbind(lphib, lphi))
      Fr <- -lc1 %*% lbphi
      lc1 <- lc1 + (1/E) * Fr %*% t(Fr)
      lc2 <- (1/E) * Fr
      lc3 <- t(lc2)
      lc4 <- matrix(1/E, 1, 1)
      lc <- cbind(rbind(lc1, lc3), rbind(lc2, lc4))
      #GLbeta[i,] <- diag(GLbetai)
      #GLphi[i,] <- diag(GLphii)
      #GLphir[i,] <- diag(GLphiri)
      #G[i,] <- diag(G)
      #GL[i,] <- GLbeta[i,] + G[i,] + GLphi[i,] + GLphir[i,]
      #Bi[i,] <- (a/dispersion[i]) * GL[i,]

  
      #DGbetai <- DGphii <- matrix(0,n,1)
      #for(j in 1:n){
      #vu2 <- V[!wzero]^2
      #Delta <- diag(0,n) ; Delta[j,j] = 1
      #Delta <- diag(1,n) - Delta; 
      #DGbetai[j] <- (vu2[j]/scale*(1 - Hi[j,j])^2)* u[j]*wi[j,j]*Hi[j,j]
      #DGphii[j] <- (sum(diag(wi))*(sqrt(vu2[j])*u[j]-1)^2*wi[j,j]^2)/((sum(diag(wi)%*%Delta))^2*scaledispersion)
      #}
      #DGbeta[i,] <- DGbetai
      #DGphi[i,] <- DGphii

      deltab <- (1/dispersion[i]) * diag((object$lm$y[!wzero] - 
                object$flm[i,]) * V[!wzero]) %*% wi %*% Xd
      deltad <- wi %*% matrix(-(0.5/dispersion[i]) * (1 - V[!wzero] * 
                u), n, 1)
      delta <- t(cbind(deltab, deltad))
      b11 <- cbind(matrix(0, p, p), matrix(0, p, 1))
      b12 <- cbind(matrix(0, 1, p), 1/E)
      b1 <- rbind(b11, b12)
      b211 <- cbind(lbb1, matrix(0, p, 1))
      b212 <- cbind(matrix(0, 1, p), matrix(0, 1, 1))
      b2 <- rbind(b211, b212)
      Cici <- -t(delta) %*% (lc) %*% delta
      Cic[[k]][i,] <- 2 * diag(Cici)
      A <- as.matrix(t(delta) %*% (lc) %*% delta)
      decA <- eigen(A)
      Lmax <- decA$val[1]
      dmax <- decA$vec[, 1]
      dmax <- dmax/sqrt(Lmax)
      dmaxc <- abs(dmax)
  
      deltab <- (-2/dispersion[i]) * t(Xd) %*% wi %*% db
      deltad <- (-1/(dispersion[i]^2)) * t(roi) %*% wi %*% db
      delta <- rbind(deltab, deltad)
      b11 <- cbind(matrix(0, p, p), matrix(0, p, 1))
      b12 <- cbind(matrix(0, 1, p), 1/lphi)
      b1 <- rbind(b11, b12)
      b211 <- cbind(lbb1, matrix(0, p, 1))
      b212 <- cbind(matrix(0, 1, p), 0)
      b2 <- rbind(b211, b212)
      Cihi <- -t(delta) %*% (lc) %*% delta
      Cih[[k]][i,] <- 2 * diag(Cihi)
      A <- as.matrix(t(delta) %*% (lc) %*% delta)
      decA <- eigen(A)
      Lmax <- decA$val[1]
      dmax <- decA$vec[, 1]
      dmax <- dmax/sqrt(Lmax)
      dmaxh <- abs(dmax)
  
      deltab <- (1/dispersion[i]) * t(Xd) %*% wi %*% da
      deltad <- (-2/(dispersion[i]^2)) * t(b) %*% wi
      delta <- rbind(deltab, deltad)
      b11 <- cbind(matrix(0, p, p), matrix(0, p, 1))
      b12 <- cbind(matrix(0, 1, p), 1/lphi)
      b1 <- rbind(b11, b12)
      b211 <- cbind(lbb1, matrix(0, p, 1))
      b212 <- cbind(matrix(0, 1, p), 0)
      b2 <- rbind(b211, b212)
      Ciri <- -t(delta) %*% (lc) %*% delta
      Cir[[k]][i,] <- 2 * diag(Ciri)
      ds <- diag(sqrt(dispersion[i]), n)
      deltai <- (1/dispersion[i]) * t(Xd) %*% da %*% wi %*% ds
      A[, i] <- as.matrix(t(deltai) %*% solve(t(Xd) %*% da %*% 
                wi %*% Xd) %*% matrix(Xd[i, ], p, 1))
	    Lmaxr[[k]][i,] <- abs(diag(A))

	    Ff <- matrix(0, p, n)
	    Ff[1, ] <- rep(1, n)
	    De <- diag(as.vector(roi), n)
	    Dv <- diag(V[!wzero])
	    st <- sqrt(var(Xd[, 1]))
	    A[, i] <- st * (1/dispersion[i]) * t(Ff %*% De %*% 
	              wi %*% Dv - object$coef[i,k] * t(Xd) %*% wi %*% da) %*% 
	    solve(t(Xd) %*% da %*% wi %*% Xd) %*% matrix(Xd[i, ], p, 1)
	    Cmaxc[k, i] <- 2 * abs(t(A[, i]) %*% A[, i])
	    Lmaxc[[k]][i, ] <- abs(diag(A))
      #for (j in 1:p) {
      #  Ff <- matrix(0, p, n)
      #  Ff[j, ] <- rep(1, n)
      #  De <- diag(as.vector(roi), n)
      #  Dv <- diag(V[!wzero])
      #  st <- sqrt(var(Xd[, j]))
      #  A[, i] <- st * (1/dispersion[i]) * t(Ff %*% De %*% 
      #            wi %*% Dv - object$coef[i,j] * t(Xd) %*% wi %*% da) %*% 
      #  solve(t(Xd) %*% da %*% wi %*% Xd) %*% matrix(Xd[i, ], p, 1)
      #  Cmaxc[j, i] <- 2 * abs(t(A[, i]) %*% A[, i])
      #  Lmaxc[[j]][i, ] <- abs(diag(A))
      #}
    }
  }
  
  Y <- object$lm$y ; yhat <- object$fitted
  res <- (Y - yhat)/sqrt(dispersion)
  rq <- rep(0,n)
  if (dist == "Normal") {
    dist.y <- pnorm(Y, yhat, dispersion)
  }
  else if (dist == "Cauchy") {
    dist.y <- pcauchy(Y, yhat, dispersion)
  }
  else if (dist == "Student") {
    dist.y <- pt(res, family$df)
  }
  else if (dist == "Gstudent") {
    dist.y <- pgstudent(res, family$s, family$r)
  }
  else if (dist == "LogisI") {
    stop(paste("not implemented yet"))
    dist.y <- plogisI(Y, yhat, dispersion)
  }
  else if (dist == "LogisII") {
    dist.y <- plogisII(res)
  }
  else if (dist == "Glogis") {
    stop(paste("not implement yet"))
    dist.y <- NULL#pglogis(res, family$alpha, family$mp)
  }
  else if (dist == "Cnormal") {
    stop(paste("not implemented yet"))
    dist.y <- NULL#pcnormal(res, family$epsi, family$sigmap)
  }
  else if (charmatch(dist, "Powerexp", F)) {
    dist.y <- ppowerexp(res, family$k)
  }
  rq <- qnorm(dist.y, 0, 1)
  
  H <- object$GW.diagnostic$Hat
  
  list(ro = ro, rr = rr, rp = rp, rq = rq, dispersion = dispersion, Hat = H, h = diag(H), 
      Cic =  lapply(Cic,function(x){diag(x)}), Cih = lapply(Cih,function(x){diag(x)}), 
      Lmaxr = lapply(Lmaxr,function(x){diag(x)}), Lmaxc = lapply(Lmaxc,function(x){diag(x)}))
}