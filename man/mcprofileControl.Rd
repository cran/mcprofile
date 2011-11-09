\name{mcprofileControl}
\alias{mcprofileControl}
\title{mcpcalc Control Arguments}
\description{Control arguments for the mcpcalc function}
\usage{
mcprofileControl(steps = 8, alphamax = 0.001, maxsteps = 200)
}
\arguments{
  \item{steps}{Number of points used for profile
    interpolation. Especially useful when defining a margin.}
  \item{alphamax}{Type-I-error used to define the maximum z-range.}
  \item{maxsteps}{Maximum number of optimization steps.}
}
\author{Daniel Gerhard}

\seealso{\code{\link{mcprofile}}}

\keyword{ iteration }
\keyword{ misc }
