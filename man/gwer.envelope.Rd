% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwer.envelope.R
\name{gwer.envelope}
\alias{gwer.envelope}
\title{Simulated Envelope of Residuals for Geographically Weighted Elliptical Regression Models}
\usage{
gwer.envelope(
  object,
  B = 100,
  arg,
  xlab = NULL,
  ylab = NULL,
  ident = NULL,
  ident.labels = NULL,
  ...
)
}
\arguments{
\item{object}{an object with the result of the fitted geographically weighted elliptical regression Models.}

\item{B}{number of Monte Carlo simulations.}

\item{arg}{a numerical or vector representing the distribution parameters used in fitted model.}

\item{xlab}{a label for the x axis, defaults to a description of x.}

\item{ylab}{a label for the y axis, defaults to a description of y.}

\item{ident}{a numerical indicate the number of observation identified in plot.}

\item{ident.labels}{an optional character vector giving labels for the identified points.}

\item{...}{arguments to be used to form the default control argument if it is not supplied directly.}
}
\description{
This function produces quantile-quantile residuals plot with simulated envelope for the specified error 
distribution in geographically weighted elliptical regression Models.
}
\examples{
\donttest{
data(georgia, package = "spgwr")
fit.formula <- PctBach ~ TotPop90 + PctRural + PctFB + PctPov
gwer.bw.t <- bw.gwer(fit.formula, data = gSRDF, family = Student(3), adapt = TRUE)
gwer.fit.t <- gwer(fit.formula, data = gSRDF, family = Student(3), bandwidth = gwer.bw.t, 
                   adapt = TRUE, parplot = FALSE, hatmatrix = TRUE, spdisp = TRUE, 
                   method = "gwer.fit")
}
}
\references{
Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
\doi{10.1016/j.spl.2005.04.033}
}
\seealso{
\code{\link{glm}}, \code{\link{elliptical}}, \code{\link{family.elliptical}}
}
\keyword{Elliptical}
\keyword{Geographically}
\keyword{Regression}
\keyword{Weighted}
\keyword{models}
\keyword{plots}
\keyword{quantile-quantile}
\keyword{regression}
