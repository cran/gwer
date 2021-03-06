% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residuals.gwer.R
\name{residuals.gwer}
\alias{residuals.gwer}
\title{Extract Residuals for Geographically Weighted Elliptical Regression Model Fits}
\usage{
\method{residuals}{gwer}(
  object,
  type = c("stand", "ordinal", "response", "pearson", "desvio"),
  ...
)
}
\arguments{
\item{object}{an object with the result of the fitted geographically weighted elliptical regression model.}

\item{type}{a character string that indicates the type of residuals. If is \code{stand} will be computed the standar residuals. 
If is \code{ordinal} will be computed the ordinal residuals. If is \code{response} will be computed the response residuals. 
If is \code{pearson} will be computed the pearson residuals. If is \code{desvio} will be computed the desviance residuals.}

\item{...}{arguments to be used to form the default control argument if it is not supplied directly.}
}
\value{
Residuals of the specific \code{type} extracted from the \code{object}.
}
\description{
This function compute different types of residuals to the fitted geographically weighted elliptical regression model.
}
\examples{
\donttest{
data(georgia, package = "spgwr")
fit.formula <- PctBach ~ TotPop90 + PctRural + PctFB + PctPov
gwer.bw.t <- bw.gwer(fit.formula, data = gSRDF, family = Student(3), adapt = TRUE)
gwer.fit.t <- gwer(fit.formula, data = gSRDF, family = Student(3), bandwidth = gwer.bw.t, 
                   adapt = TRUE, parplot = FALSE, hatmatrix = TRUE, spdisp = TRUE, 
                   method = "gwer.fit")
summary(gwer.fit.t) 
residuals(gwer.fit.t, type = "stand") 
}
}
\references{
Brunsdon, C., Fotheringham, A. S. and Charlton, M. E. (1996). 
Geographically weighted regression: a method for exploring spatial nonstationarity.
Geographical analysis, 28(4), 281-298. \doi{10.1111/j.1538-4632.1996.tb00936.x}

Galea, M., Paula, G. A., and Cysneiros, F. J. A. (2005). On diagnostics in 
symmetrical nonlinear models. Statistics & Probability Letters, 73(4), 459-467.
\doi{10.1016/j.spl.2005.04.033}
}
\seealso{
\code{\link{residuals}}, \code{\link{gwer}}, \code{\link{family.elliptical}}
}
\keyword{Elliptical}
\keyword{Geographically}
\keyword{Residuals}
\keyword{models}
\keyword{regression}
\keyword{weighted}
