\name{confint-methods}
\docType{methods}
\alias{confint,mcprofile-method}
\alias{confint,mcprofileRatio-method}
\title{Simultaneous Confidence Intervals for Multiple Contrast Profiles}
\description{
 Calculates simultaneous confidence intervals based on
 \code{\linkS4class{mcprofile}} objects. Confidence limits are found by
   spline interpolation at a signed root deviance cutoff value.
   Multiplicity adjustment can be performed by Bonferroni correction or
   utilizing a multivariate normal or t-distribution.
}
\usage{
\S4method{confint}{mcprofile}(object, parm, level=0.95,
adjust="single-step", alternative="two.sided", quant=NULL)
}

\arguments{
  \item{object}{An object of class \code{\linkS4class{mcprofile}}, or \code{\linkS4class{mcprofileRatio}}}
  \item{parm}{Missing, only needed to provide a mcprofile method for \code{confint}}
  \item{level}{Global confidence level, default at 0.95}
  \item{adjust}{Multiplicity adjustment method: \code{"none"} for unadjusted
  intervals, \code{"bonferroni"}, or \code{"single-step"} for utilizing a multivariate
  normal or t-distribution.}
  \item{alternative}{Choosing one- or two-sided intervals by
  \code{"two.sided"}, \code{"less"}, or \code{"greater"}}
  \item{quant}{Numeric, in not \code{NULL} this value is used as cutoff
  value for the signed root deviance}
}

\details{Confidence limits are found by spline interpolation using the
  splines in the \code{bsplines} slot of the
  \code{\linkS4class{mcprofile}}, or \code{\linkS4class{mcprofileRatio}} object.  Cutoff values are calculated
  as quantiles of a normal distribution, or t-distribution if the df
  slot is not empty (gaussian families).  For \code{"single-step"}
  adjustment the correlation structure of a multivariate normal
  distribution is estimated by standardizing the variance-covariance
  matrix of the linear combinations of parameters.  This vcov matrix is
  obtained by multiplying the original vcov matrix with the prespecified
  contrast matrix from both sides.}

\value{
An object of class \code{\linkS4class{mcpconfint}}
}

\seealso{\code{\link{mcpcalc}}, \code{\link{mcpcalcRatio}}, \code{\link{confint.glht}}, \code{\link{sci.ratio}}}

\keyword{methods}
\keyword{htest}
