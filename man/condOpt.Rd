\name{condOpt}
\alias{condOpt}
\title{Conditional Optimization}
\description{
  Performs conditional optimization given an object constructed by
  function glmobj, a designmatrix and an offset. Only for internal use.
}
\usage{
condOpt(object, X, offset, start=NULL)
}
\arguments{
  \item{object}{an object constructed by function glmobj.}
  \item{X}{a design matrix}
  \item{offset}{a numeric offset vector}
  \item{start}{a numeric vector with starting values}
}
\details{
  Uses \code{\link{glm.fit}} to perform conditional optimization.
}
\keyword{optimize}
