#' @title Diagnostic for Geographically Weighted Elliptical Regression Models
#' @description This function obtains the values of differents residuals types and calculates the diagnostic measures for the fitted geographically weighted elliptical regression model.
#' @param object an object with the result of the fitted geographically weighted elliptical regression model.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return Returns a list of diagnostic arrays:
#' \item{ro}{Ordinal residuals.}
#' \item{rr}{Response residuals.}
#' \item{rp}{Pearson residuals.}
#' \item{rs}{Studentized residuals.}
#' \item{rd}{Deviance residuals.}
#' \item{dispersion}{Coefficient of dispersion.}
#' \item{H}{The hat matrix.}
#' \item{h}{Main diagonal of the hat matrix.}
#' \item{GL}{Generalized leverage.}
#' \item{GLbeta}{Generalized leverage of location parameters estimation.}
#' \item{GLphi}{Generalized leverage of dispersion parameters estimation.}  
#' \item{DGbeta}{Cook distance of location parameters estimation.}
#' \item{DGphi}{Cook distance of dispersion parameters estimation.}  
#' \item{Cic}{Normal curvature for case-weight perturbation.}
#' \item{Cih}{Normal curvature for scale perturbation.}
#' \item{Lmaxr}{Local influence on response (additive perturbation in responce).}
#' \item{Lmaxc}{Local influence on coefficients (additive perturbation in predictors).}
#' @references Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
#' symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
#' \url{https://doi.org/10.1016/j.spl.2005.04.033}
#' @seealso \code{\link{elliptical}}
#' @keywords Geographically Weighted Regression
#' @keywords Elliptical models
#' @keywords Diagnostic methods
#' @examples
#' data(columbus, package="spData")
#' fit.lm <- lm(CRIME ~ INC, data=columbus)
#' summary(fit.lm)
#' gwer.bw <- gwer.sel(CRIME ~ INC, data=columbus, family = Normal(),
#'                  coords=cbind(columbus$X, columbus$Y))
#' gwer.fitn <- gwer(CRIME ~ INC, family = Normal(), bandwidth = gwer.bw, hatmatrix = TRUE, 
#'                  spdisp = TRUE, parplot = TRUE, data=columbus, method = "gwer.fit",
#'                  coords=cbind(columbus$X, columbus$Y))
#' gwer.diag(gwer.fitn)               
#' \donttest{
#' data(columbus, package="spData")
#' fit.elliptical <- elliptical(CRIME ~ INC, family = Student(df=3), data=columbus)
#' summary(fit.elliptical)
#' gwer.bw <- gwer.sel(CRIME ~ INC, data=columbus, family = Student(df=3),
#'                  coords=cbind(columbus$X, columbus$Y), method = 'aic')
#' gwer.fitt <- gwer(CRIME ~ INC, family = Student(df=3), bandwidth = gwer.bw, 
#'                  spdisp = TRUE, hatmatrix = TRUE, data=columbus, method = "gwer.fit",
#'                  coords=cbind(columbus$X, columbus$Y))
#' gwer.diag(gwer.fitt)  
#' }
#' @export

gwer.diag <- function (object,...) 
{
  if (object$fp.given || !object$hatmatrix) 
    stop("Diagnostic measures not applicable - regression points different from observed points")
  
  Wi <- object$gweights ; l.fitted = object$flm
  p <- object$lm$rank
  family <- object$family
  user.def <- object$lm$user.def
  f.name <- family[[1]]
  dispersion <- as.numeric(object$dispersion)
  w <- if (is.null(object$lm$weights)) 
    rep(1, length(object$lm$residuals))
  else object$lm$weights
  wzero <- (w == 0)
  Xd <- diag(c(w[!wzero])) %*% object$lm$Xmodel[!wzero, 
                                                       ]
  n <- length(object$lm$residuals)
  if(length(dispersion)==1)
    dispersion <- rep(dispersion,n)
  ro <- rr <- rp <- rs <- rd <- rep(0,n)
  H <- G <- GL <- GLbeta <- GLphi <- GLphir <- Bi <- Lmaxr <- Cic <- Cih <- DGbeta <- DGphi <- DGphir <- matrix(0,n,n)
  Cmaxc <- Lmaxc <- matrix(0, p, n) ; Lmaxc <- list() ; for(j in 1:p){Lmaxc[[j]] <- matrix(0,n,n)}
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
    Hi <- Xd %*% solve(t(Xd) %*% wi %*% Xd) %*% t(Xd) %*% wi
    H[i,] <- Xd[i,] %*% solve(t(Xd) %*% wi %*% Xd) %*% t(Xd) %*% wi
    h <- diag(Hi)/(scalevariance * scale)
    rs[i] <- resid[i]/sqrt(scalevariance * (1 - h[i]))

    logg0 <- family$g0(0, df = family$df, s = family$s, 
                       r = family$r, alpha = family$alpha, mp = family$mp, 
                       epsi = family$epsi, sigmap = family$sigmap, k = family$k)
    loggu <- family$g0(u[i], df = family$df, s = family$s, 
                       r = family$r, alpha = family$alpha, mp = family$mp, 
                       epsi = family$epsi, sigmap = family$sigmap, k = family$k)
    rd[i] <- sqrt(Wi[i,i])*(sign(resid[i]))*(2*logg0 - 2*loggu)^(.5)
  

    ct <- Wgder[!wzero]
    op <- Wg[!wzero] * roi
    a <- V[!wzero] - 4 * ct * u
    b <- op + (u * ct * roi)
    db <- diag(b)
    b <- matrix(b, n, 1)
    u <- matrix(u, n, 1)
    da <- diag(a)
    dai <- diag(1/a)
    roi <- matrix(roi, n, 1)
    som <- matrix(0, p, p)
    Gi <- GLphiri <- matrix(0, n, n)
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
    Fi <- -lbb1 %*% lbphi
    GLbetai <- Xd %*% (-lbb1) %*% lby
    R <- Xd %*% (Fi %*% t(Fi) %*% lby + Fi %*% lphiy)
    GLphii <- (-1/E) * R
  
    Om <- rbind(cbind(lbb, lbphi), cbind(lphib, lphi))
    Fr <- -lc1 %*% lbphi
    lc1 <- lc1 + (1/E) * Fr %*% t(Fr)
    lc2 <- (1/E) * Fr
    lc3 <- t(lc2)
    lc4 <- matrix(1/E, 1, 1)
    lc <- cbind(rbind(lc1, lc3), rbind(lc2, lc4))
    GLbeta[i,] <- diag(GLbetai)
    GLphi[i,] <- diag(GLphii)
    GLphir[i,] <- diag(GLphiri)
    G[i,] <- diag(G)
    GL[i,] <- GLbeta[i,] + G[i,] + GLphi[i,] + GLphir[i,]
    Bi[i,] <- (a/dispersion[i]) * GL[i,]

  
    DGbetai <- DGphii <- matrix(0,n,1)
    for(j in 1:n){
    vu2 <- V[!wzero]^2
    Delta <- diag(0,n) ; Delta[j,j] = 1
    Delta <- diag(1,n) - Delta; 
    DGbetai[j] <- (vu2[j]/scale*(1 - Hi[j,j])^2)* u[j]*wi[j,j]*Hi[j,j]
    DGphii[j] <- (sum(diag(wi))*(sqrt(vu2[j])*u[j]-1)^2*wi[j,j]^2)/((sum(diag(wi)%*%Delta))^2*scaledispersion)
    }
    DGbeta[i,] <- DGbetai
    DGphi[i,] <- DGphii

    deltab <- (1/dispersion[i]) * diag((object$lm$y[!wzero] - 
                                     object$flm[i,]) * V[!wzero]) %*% 
      wi %*% Xd
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
    Cic[i,] <- 2 * diag(Cici)
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
    Cih[i,] <- 2 * diag(Cihi)
    A <- as.matrix(t(delta) %*% (lc) %*% delta)
    decA <- eigen(A)
    Lmax <- decA$val[1]
    dmax <- decA$vec[, 1]
    dmax <- dmax/sqrt(Lmax)
    dmax <- abs(dmax)
  
    deltab <- (1/dispersion[i]) * t(Xd) %*% wi %*% da
    deltad <- (-2/(dispersion[i]^2)) * t(b) %*% wi
    delta <- rbind(deltab, deltad)
    b11 <- cbind(matrix(0, p, p), matrix(0, p, 1))
    b12 <- cbind(matrix(0, 1, p), 1/lphi)
    b1 <- rbind(b11, b12)
    b211 <- cbind(lbb1, matrix(0, p, 1))
    b212 <- cbind(matrix(0, 1, p), 0)
    b2 <- rbind(b211, b212)
    Ci <- -t(delta) %*% (lc) %*% delta
    Ci <- 2 * diag(Ci)

    ds <- diag(sqrt(dispersion[i]), n)
    deltai <- (1/dispersion[i]) * t(Xd) %*% da %*% wi %*% ds
    A[, i] <- as.matrix(t(deltai) %*% solve(t(Xd) %*% da %*% 
                                              wi %*% Xd) %*% matrix(Xd[i, ], p, 1))
	
    Lmaxr[i,] <- abs(diag(A))

    for (j in 1:p) {
      Ff <- matrix(0, p, n)
      Ff[j, ] <- rep(1, n)
      De <- diag(as.vector(roi), n)
      Dv <- diag(V[!wzero])
      st <- sqrt(var(Xd[, j]))
      A[, i] <- st * (1/dispersion[i]) * t(Ff %*% De %*% 
                                          wi %*% Dv - object$coef$est[i,j] * t(Xd) %*% wi %*% da) %*% 
        solve(t(Xd) %*% da %*% wi %*% Xd) %*% matrix(Xd[i, 
                                                        ], p, 1)
      Cmaxc[j, i] <- 2 * abs(t(A[, i]) %*% A[, i])
      Lmaxc[[j]][i, ] <- abs(diag(A))
    }
  }
  
  
  list(ro = ro, rr = rr, rp = rp, rs = rs, rd = rd, dispersion = dispersion, 
       Hat = H, h = diag(H), DGbeta = diag(DGbeta), DGphi = diag(DGphi), 
       GL = diag(GL), GLbeta = diag(GLbeta), GLphi = diag(GLphi), 
       Cic = diag(Cic), Cih = diag(Cih), Lmaxr = diag(Lmaxr), 
       Lmaxc = lapply(Lmaxc,function(x){diag(x)}))
}