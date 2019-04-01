#' @title Diagnostic for Elliptical Regression Models
#' @description This function obtains the values of the residuals and calculates the diagnostic measures for elliptical regression models.
#' @param ellipticalfit fit object for elliptical regression model.
#' @param weighting type of model weighting used.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#' @return returns a list of diagnostic arrays:
#' \item{resid}{ordinal residuals for the fit model.}
#' \item{rs}{studentized residuals for the fit model.}
#' \item{dispersion}{coefficient of dispersion for the model fit.}
#' \item{GL}{generalized leverage for the model fit.}
#' \item{GLbeta}{generalized leverage of location parameters estimation for the model fit.}
#' \item{GLphi}{generalized leverage of dispersion parameters estimation for the model fit.}  
#' \item{Bi}{generalized leverage weighted by dispersion for the model fit.}
#' \item{Om}{observed fisher information matrix of the model fit.}
#' \item{Iom}{expected fisher information matrix of the model fit.}  
#' \item{a, b, c}{the value of D(a), D(b) and D(c), respectively, for the model fit.}
#' \item{Cmax}{matrix of local influence for additive perturbation in response.}
#' \item{Lmax}{matrix of local influence on coefficients (additive perturbation in predictors).}
#' \item{Cic}{matrix of local influence for case-weight perturbation (Ci).}
#' \item{dmax}{matrix of local influence for case-weight perturbation (dmax).}
#' \item{dmaxc}{matrix of local influence for case-weight perturbation (|dmax|).}
#' \item{Ci}{matrix of local influence on the scale.}
#' \item{Cih}{main diagonal of the matrix of local influence on the scale.}
#' \item{h}{main diagonal of the hat matrix.}
#' @references Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
#' symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
#' \url{https://doi.org/10.1016/j.spl.2005.04.033}
#' @seealso \code{\link{elliptical}}
#' @keywords elliptical
#' @keywords diagnostic methods
#' @examples
#'data(luzdat)
#' y <- luzdat$y
#' x1 <- luzdat$x1 ; x1 <- factor(x1) ; x1 <- C(x1,treatment)
#' x2 <- luzdat$x2
#' x3 <- (luzdat$x2)^2
#' luz <- data.frame(y,x1,x2,x3)
#' elliptical.fitt <- elliptical(y ~ x1+x2+x3, family = Student(df=5),
#' data=luz)
#' elliptical.diag(elliptical.fitt)
#' @export

elliptical.diag <- function (ellipticalfit, weighting = "observed",...) 
{
  scalevariance <- ellipticalfit$scalevariance
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
  H <- Xd %*% solve(t(Xd) %*% Xd) %*% t(Xd)
  h <- diag(H)/(scalevariance * scale)
  rs <- resid/sqrt(scalevariance * (1 - h))
  ro <- ellipticalfit$y[!wzero] - ellipticalfit$fitted[!wzero]
  n <- length(resid)
  u <- ro^2/dispersion
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
  M <- as.matrix(som + t(Xd) %*% da %*% Xd)
  lbb <- (-1/dispersion) * M
  lby <- (1/dispersion) * t(Xd) %*% da
  lphiy <- (-2/(dispersion^2)) * t(b)
  lbphi <- (2/dispersion^2) * t(Xd) %*% b
  lphib <- t(lbphi)
  lphi <- (1/(dispersion^2)) * ((n/2) + t(u) %*% diag(ct) %*% 
                                  u - (1/dispersion) * t(ro) %*% diag(ellipticalfit$v[!wzero]) %*% 
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
  deltab <- matrix(0, n, p)
  deltad <- matrix(0, n, 1)
  deltab <- (1/dispersion) * diag((ellipticalfit$y[!wzero] - 
                                     ellipticalfit$fitted[!wzero]) * ellipticalfit$v[!wzero]) %*% 
    Xd
  deltad <- matrix(-(0.5/dispersion) * (1 - ellipticalfit$v[!wzero] * 
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
  Lmax <- NULL
  A <- matrix(0, n, n)
  for (i in 1:n) {
    A[, i] <- as.matrix(t(deltai) %*% solve(t(Xd) %*% da %*% 
                                              Xd) %*% matrix(Xd[i, ], p, 1))
  }
  Lmax <- abs(diag(A))
  Cmax <- matrix(0, p, n)
  for (j in 1:p) {
    Ff <- matrix(0, p, n)
    Ff[j, ] <- rep(1, n)
    De <- diag(as.vector(ro), n)
    Dv <- diag(ellipticalfit$v[!wzero])
    st <- sqrt(var(Xd[, j]))
    for (i in 1:n) {
      A[, i] <- st * (1/dispersion) * t(Ff %*% De %*% 
                                          Dv - ellipticalfit$coef[j] * t(Xd) %*% da) %*% 
        solve(t(Xd) %*% da %*% Xd) %*% matrix(Xd[i, 
                                                 ], p, 1)
      Cmax[j, i] <- 2 * abs(t(A[, i]) %*% A[, i])
    }
  }
  list(resid = resid, rs = rs, dispersion = dispersion, GL = GL, 
       GLbeta = GLbeta, GLphi = GLphi, G = G, GLphir = GLphir, 
       dmax = dmax, Ci = Ci, Bi = Bi, Om = Om,  IOm = lc, a = a,
       b = as.vector(b), c = ct, Cmax = Cmax, Lmax = Lmax,
       Cic = Cic, dmaxc = dmaxc, Cih = Cih, h = diag(H))
}
