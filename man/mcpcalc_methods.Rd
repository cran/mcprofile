\name{mcpcalc_methods}
\alias{mcpcalcBFGS}
\alias{mcpcalcIRWLS}
\alias{mcpcalcnls}
\alias{mcpcalcmer}
\alias{mcpcalcmerlmer}

\title{Construction of Multiple Contrast Profiles} 
\description{Only for internal use. Calculating profiles by rearranging
  the design matrix (IRWLS, lmer) or with an optim wrapper using BFGS method}

\usage{
mcpcalcIRWLS(object, CM, control=mcprofileControl(), margin=NULL)
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

\details{
  According to the rows of a contrast matrix the design matrix of
a model is restructured by \code{\link{makeDesign}} (IRWLS, lmer). For
  several constrained model fits the signed root deviance is calculated and interpolated by interpolation splines.  
}

\value{ An object inheriting from class \code{\linkS4class{mcprofile}}
  with methods \code{\link{show}}, \code{\link{plot}},
  \code{\link{confint}}, and \code{\link{test}}}

\seealso{\code{\link{profile.glm}}, \code{\link{glht}}, \code{\link{interpSpline}}, \code{\link{contrMat}}}

\keyword{models}
\keyword{optimize}

