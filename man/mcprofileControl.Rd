\name{mcprofileControl}
\alias{mcprofileControl}
\title{Constructs a list with mcpcalc control parameters.} 
\description{Constructs a list for the control argument in a mcpcalc
  function call.}

\usage{
mcprofileControl(steps=8, alphamax=0.001, fixed.range=FALSE, maxsteps=200)
}

\arguments{
  \item{steps}{Step size for profile construction, or number of deviance
    evaluations if fixed.range=TRUE.}
  \item{alphamax}{Maximal alpha level to define maximum profile
    statistics by choosing a quantile of a normal distribution at 1-alphamax 
divided by the number of rows of the contrast matrix.}
  \item{fixed.range}{If TRUE, the profile is calculated for equidistant
    points between two margins.}
  \item{maxsteps}{Maximum number of points for profile calculation.}
}

\value{
  A list control parameters.
}
\keyword{models}
