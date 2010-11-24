\name{test-methods}
\docType{methods}
\alias{test}
\alias{test-methods}
\alias{test,mcprofile-method}
\alias{test,mcprofileRatio-method}

\title{Multiple Tests for Multiple Contrast Profiles}
\description{
  Calculates adjusted p-values based on \code{\linkS4class{mcprofile}},
  or \code{\linkS4class{mcprofileRatio}} objects. Test statistics are found by spline interpolation at a specified test margin. Multiplicity adjustment can be performed by utilizing a multivariate normal or t-distribution, or by any method
  provided by \code{\link{p.adjust}}.
}

\usage{
\S4method{test}{mcprofile}(object, adjust="single-step",
alternative="two.sided", margin=0)
}

\arguments{
  \item{object}{An object of class \code{\linkS4class{mcprofile}}, or \code{\linkS4class{mcprofileRatio}}}
  \item{adjust}{Multiplicity adjustment method: \code{"single-step"} for
  utilizing a multivariate normal or t-distribution, or any method
  provided by \code{\link{p.adjust}}.}
  \item{alternative}{One- or two-sided testing by \code{"two.sided"}, \code{"less"}, or \code{"greater"}}
  \item{margin}{Test margin to define a null hypothesis; until now only
  a single margin for all tests. If the profile range does not include
  the test margin, test statistics are found by extrapolation!}
}

\details{Test statistics are found by spline interpolation using the
  splines in the \code{fsplines} slot of the
  \code{\linkS4class{mcprofile}}, or \code{\linkS4class{mcprofileRatio}} object.  p-values are calculated from a
  normal distribution, or t-distribution if the \code{df} slot is not empty (gaussian families).  For \code{"single-step"} adjustment the correlation structure of a multivariate normal distribution is estimated by standardizing the variance-covariance matrix of the linear combinations of parameters.  This vcov matrix is obtained by multiplying the original vcov matrix with the prespecified contrast matrix from both sides.}

\value{
An object of class \code{\linkS4class{mcptest}}
}

\seealso{\code{\link{mcpcalc}}, \code{\link{glht}}, \code{\link{simtest.ratio}}}


\keyword{methods}
\keyword{htest}
