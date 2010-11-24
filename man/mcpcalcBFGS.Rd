\name{mcpcalcBFGS}
\alias{mcpcalcBFGS}
\title{Construction of Multiple Contrast Profiles} 
\description{Only for internal use. Calculating profiles by optimization
with linear constraints using optim.}

\usage{
mcpcalcBFGS(object, CM, control=mcprofileControl(), margin=NULL)
}

\arguments{
  \item{object}{An object of class \code{\link{glm}} or \code{\link{lm}}}
  \item{CM}{A contrast matrix with the number of columns equal the
  number of estimated parameters.  Contrast coefficients should be -1
  and 1, respectively.  Providing row names is recommendable.}
  \item{control}{A list with control arguments. See \code{\link{mcprofileControl}}.}
  \item{margin}{Numeric vector of margins, defining the profiling
  range.  If \code{NULL} no extra margin is considered.}
}

\details{ The signed root deviance is calculated and interpolated by interpolation splines.  Control arguments are provided to define the profiling range:
\code{steps} correspond to the number of values evaluated around the
MLE, \code{alphamax} defines the maximal cutoff value by choosing a
quantile of a normal distribution at 1-alphamax divided by the number of
rows of the contrast matrix, if \code{fixed.range} is \code{TRUE} the
deviance is evaluated at equal-sized steps in a range provided by the
\code{margin} argument.
}

\value{ An object inheriting from class \code{\linkS4class{mcprofile}}
  with methods \code{\link{show}}, \code{\link{plot}},
  \code{\link{confint}}, and \code{\link{test}}}

\seealso{\code{\link{optim}}, \code{\link{profile.glm}}, \code{\link{glht}}, \code{\link{interpSpline}}, \code{\link{contrMat}}}

\keyword{models}
\keyword{optimize}

