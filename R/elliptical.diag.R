#' @title Diagnostic Measures for Elliptical Regression Models
#' @description This function obtains the values of different residuals types and calculates the diagnostic measures for the fitted elliptical regression model.
#' @param object an object with the result of the fitted elliptical regression model.
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
#' @references Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
#' symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
#' \doi{10.1016/j.spl.2005.04.033}
#' @seealso \code{\link{elliptical}}
#' @keywords Elliptical regression models
#' @keywords Diagnostic methods
#' @examples
#' data(luzdat)
#' y <- luzdat$y
#' x1 <- luzdat$x1 ; x1 <- factor(x1) ; x1 <- C(x1, treatment)
#' x2 <- luzdat$x2
#' x3 <- (luzdat$x2)^2
#' luz <- data.frame(y, x1, x2, x3)
#' elliptical.fitt <- elliptical(y ~ x1+x2+x3, family = Student(df=5),
#' data = luz)
#' elliptical.diag(elliptical.fitt)
#' @export

elliptical.diag <- function (object,...) 
{
  scalevariance <- object$scalevariance
  scaledispersion <- object$scaledispersion
  scale <- object$scale
  family <- object$family
  user.def <- object$user.def
  f.name <- family[[1]]
  dispersion <- object$dispersion
  w <- if (is.null(object$weights)) 
    rep(1, length(object$residuals))
  else object$weights
  wzero <- (w == 0)
  Xd <- diag(c(w[!wzero])) %*% object$Xmodel[!wzero, 
                                                    ]
  n <- length(object$residuals)
  ro <- object$y[!wzero] - object$fitted[!wzero]
  rr <- (object$y[!wzero] - object$fitted[!wzero])/sqrt(dispersion)
  Ones <- rep(1,n) ; u <- ro^2/dispersion
  rp <- rr/sqrt(scalevariance)

  dev <- 2 * (family$g0(rr, df = family$df, r = family$r, 
                        s = family$s, alpha = family$alpha, mp = family$mp, 
                        epsi = family$epsi, sigmap = family$sigmap, k = family$k) - 
                family$g0(0, df = family$df, r = family$r, s = family$s, 
                          alpha = family$alpha, mp = family$mp, epsi = family$epsi, 
                          sigmap = family$sigmap, k = family$k))
  p <- object$rank
  H <- Xd %*% solve(t(Xd) %*% Xd) %*% t(Xd)
  h <- diag(H)/(scalevariance * scale)
  rs <- rr/sqrt(scalevariance * (1 - h))

  
  rd <- rep(0,n)
  for(i in 1:n){
    logg0 <- family$g0(0, df = family$df, s = family$s, 
                       r = family$r, alpha = family$alpha, mp = family$mp, 
                       epsi = family$epsi, sigmap = family$sigmap, k = family$k)
    loggu <- family$g0(u[i], df = family$df, s = family$s, 
                       r = family$r, alpha = family$alpha, mp = family$mp, 
                       epsi = family$epsi, sigmap = family$sigmap, k = family$k)
    rd[i] <- (sign(rr[i]))*(2*logg0 - 2*loggu)^(.5)
  }
  
  ct <- object$Wgder[!wzero]
  op <- object$Wg[!wzero] * ro
  a <- object$v[!wzero] - 4 * ct * u
  b <- op + (u * ct * ro)
  db <- diag(b)
  b <- matrix(b, n, 1)
  u <- matrix(u, n, 1)
  da <- diag(a)
  dai <- diag(1/a)
  ro <- matrix(ro, n, 1)
  som <- matrix(0, p, p)
  M <- as.matrix(som + t(Xd) %*% da %*% Xd)
  lbb <- (-1/dispersion) * M
  lby <- (1/dispersion) * t(Xd) %*% da
  lphiy <- (-2/(dispersion^2)) * t(b)
  lbphi <- (2/dispersion^2) * t(Xd) %*% b
  lphib <- t(lbphi)
  lphi <- (1/(dispersion^2)) * ((n/2) + t(u) %*% diag(ct) %*% 
                                  u - (1/dispersion) * t(ro) %*% diag(object$v[!wzero]) %*% 
                                  ro)
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
  vu2 <- object$v[!wzero]^2
  for(i in 1:n){
    Delta <- diag(0,n) ; Delta[i,i] = 1
    Delta <- diag(1,n) - Delta; 
    DGbeta[i] <- (vu2[i]/scale*(1 - H[i,i])^2)* u[i]*H[i,i]
    DGphi[i] <- n*(sqrt(vu2[i])*u[i]-1)^2/((n-1)^2*scaledispersion)
  }
  
  deltab <- matrix(0, n, p)
  deltad <- matrix(0, n, 1)
  deltab <- (1/dispersion) * diag((object$y[!wzero] - 
                                     object$fitted[!wzero]) * object$v[!wzero]) %*% 
    Xd
  deltad <- matrix(-(0.5/dispersion) * (1 - object$v[!wzero] * 
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
  
  deltab <- (-2/dispersion) * t(Xd) %*% db
  deltad <- (-1/(dispersion^2)) * t(ro) %*% db
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
  
  deltab <- (1/dispersion) * t(Xd) %*% da
  deltad <- (-2/(dispersion^2)) * t(b)
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
  deltai <- (1/dispersion) * t(Xd) %*% da %*% ds
  Lmaxr <- NULL
  A <- matrix(0, n, n)
  for (i in 1:n) {
    A[, i] <- as.matrix(t(deltai) %*% solve(t(Xd) %*% da %*% 
                                              Xd) %*% matrix(Xd[i, ], p, 1))
  }
  Lmaxr <- abs(diag(A))
  
  Cmaxc <- Lmaxc <- matrix(0, p, n)
  for (j in 1:p) {
    Ff <- matrix(0, p, n)
    Ff[j, ] <- rep(1, n)
    De <- diag(as.vector(ro), n)
    Dv <- diag(object$v[!wzero])
    st <- sqrt(var(Xd[, j]))
    for (i in 1:n) {
      A[, i] <- st * (1/dispersion) * t(Ff %*% De %*% 
                                          Dv - object$coef[j] * t(Xd) %*% da) %*% 
        solve(t(Xd) %*% da %*% Xd) %*% matrix(Xd[i, 
                                                 ], p, 1)
      Cmaxc[j, i] <- 2 * abs(t(A[, i]) %*% A[, i])
    }
    Lmaxc[j, ] <- abs(diag(A))
  }
  list(ro = ro , rr = rr, rp = rp, rs = rs, rd = rd, dispersion = dispersion, Hat = H, 
       h = diag(H), GL = GL, GLbeta = GLbeta, GLphi = GLphi, DGbeta = DGbeta, DGphi = DGphi, 
       Cic = Cic, Cih = Cih, Lmaxr = Lmaxr, Lmaxc = Lmaxc)
}
