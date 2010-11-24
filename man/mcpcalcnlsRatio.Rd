\name{mcpcalcnlsRatio}
\alias{mcpcalcnlsRatio}
\title{Construction of Multiple Contrast Profiles} 
\description{Only for internal use. Calculating profiles for nls objects
  by optimization with linear constraints using optim.}

\usage{
mcpcalcnlsRatio(object, CMn, CMd, control=mcprofileControl(), margin=NULL)
}

\arguments{
  \item{object}{An object of class \code{\link{nls}}}
  \item{CMn}{A numerator contrast matrix with the number of columns equal the
  number of estimated parameters.  Contrast coefficients should be -1
  and 1, respectively.  Providing row names is recommendable.}
  \item{CMd}{A denominator contrast matrix equivalent to the numerator matrix.}  
  \item{control}{A list with control arguments. See the details section.}
  \item{margin}{Numeric vector of margins, defining the profiling
  range.  If \code{NULL} no extra margin is considered.}
}

\details{ The signed root deviance is calculated and interpolated by interpolation splines. 
}

\value{ An object inheriting from class \code{\linkS4class{mcprofileRatio}}
  with methods \code{\link{show}}, \code{\link{plot}},
  \code{\link{confint}}, and \code{\link{test}}}

\seealso{\code{\link{optim}}, \code{\link{profile.nls}}, \code{\link{sci.ratio}}, \code{\link{interpSpline}}, \code{\link{contrMatRatio}}}

\keyword{models}
\keyword{optimize}

