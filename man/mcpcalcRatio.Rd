\name{mcpcalcRatio}
\alias{mcpcalcRatio}
\alias{mcpcalcRatio.lm}
\alias{mcpcalcRatio.glm}
\alias{mcpcalcRatio.nls}
\title{Construction of Multiple Contrast Profiles for Ratios of Parameters} 
\description{Calculates signed root deviance profiles given a \code{\link{glm}} or \code{\link{lm}} object.  The profiled parameters of interest are defined by providing a contrast matrix.  }

\usage{
\method{mcpcalcRatio}{lm}(object, CMn, CMd, control=mcprofileControl(), margin=NULL, method="BFGS")
\method{mcpcalcRatio}{glm}(object, CMn, CMd, control=mcprofileControl(), margin=NULL, method="BFGS")
\method{mcpcalcRatio}{nls}(object, CMn, CMd, control=mcprofileControl(), margin=NULL, method="BFGS")
}

\arguments{
  \item{object}{An object of class \code{\link{glm}} or \code{\link{lm}}}
  \item{CMn}{A numerator contrast matrix with the number of columns
  equal to the number of estimated parameters.  Contrast coefficients should be -1
  and 1, respectively.  Providing row names is recommendable.}
  \item{CMd}{A denominator contrast matrix equivalent to the numerator
  matrix.}
  \item{control}{A list with control arguments. See \code{\link{mcprofileControl}}.}
  \item{margin}{Numeric vector of margins, defining the profiling
  range.  If \code{NULL} no extra margin is considered.}
  \item{method}{"BFGS"}
}

\details{
  The signed root deviance is calculated for the ratio of parameters and interpolated by interpolation splines.
}

\value{ An object inheriting from class \code{\linkS4class{mcprofileRatio}}
  with methods \code{\link{show}}, \code{\link{plot}},
  \code{\link{confint}}, and \code{\link{test}}}

\seealso{\code{\link{profile.glm}}, \code{\link{profile.nls}}, \code{\link{sci.ratio}}, \code{\link{interpSpline}}, \code{\link{contrMatRatio}}}



\keyword{models}
\keyword{optimize}

