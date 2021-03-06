% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elliptical.diag.R
\name{elliptical.diag}
\alias{elliptical.diag}
\title{Diagnostic Measures for Elliptical Regression Models}
\usage{
elliptical.diag(object, ...)
}
\arguments{
\item{object}{an object with the result of the fitted elliptical regression model.}

\item{...}{arguments to be used to form the default control argument if it is not supplied directly.}
}
\value{
Returns a list of diagnostic arrays:
\item{ro}{ordinal residuals.}
\item{rr}{response residuals.}
\item{rp}{pearson residuals.}
\item{rs}{studentized residuals.}
\item{rd}{deviance residuals.}
\item{dispersion}{coefficient of dispersion parameter.}
\item{Hat}{the hat matrix.}
\item{h}{main diagonal of the hat matrix.}
\item{GL}{generalized leverage.}
\item{GLbeta}{generalized leverage of location parameters estimation.}
\item{GLphi}{generalized leverage of dispersion parameters estimation.}  
\item{DGbeta}{cook distance of location parameters estimation.}
\item{DGphi}{cook distance of dispersion parameters estimation.}  
\item{Cic}{normal curvature for case-weight perturbation.}
\item{Cih}{normal curvature for scale perturbation.}
\item{Lmaxr}{local influence on response (additive perturbation in responce).}
\item{Lmaxc}{local influence on coefficients (additive perturbation in predictors).}
}
\description{
This function obtains the values of different residuals types and calculates the diagnostic measures for the fitted elliptical regression model.
}
\examples{
data(luzdat)
y <- luzdat$y
x1 <- luzdat$x1 ; x1 <- factor(x1) ; x1 <- C(x1, treatment)
x2 <- luzdat$x2
x3 <- (luzdat$x2)^2
luz <- data.frame(y, x1, x2, x3)
elliptical.fitt <- elliptical(y ~ x1+x2+x3, family = Student(df=5),
data = luz)
elliptical.diag(elliptical.fitt)
}
\references{
Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
\doi{10.1016/j.spl.2005.04.033}
}
\seealso{
\code{\link{elliptical}}
}
\keyword{Diagnostic}
\keyword{Elliptical}
\keyword{methods}
\keyword{models}
\keyword{regression}
