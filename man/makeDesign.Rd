\name{makeDesign}
\alias{makeDesign}
\alias{pmakeDesign}
\title{Construction of a modified design matrix. Only for internal use.} 
\description{Given a contrast and a design matrix, a new designmatrix is
constructed to enable the calculation of linear combinations of the
original model parameters.}

\usage{
makeDesign(k, X)
pmakeDesign(k, X)
}
\arguments{
  \item{k}{a numeric vector with contrast coefficients}
  \item{X}{a design matrix}
}
\details{
 If both negative and positive contrast coefficients are provided,
 \code{makeDesign} sums up all design matrix coefficients corresponding
 to non-negative contrast coefficients as a new intercept.  A dummy
 variable for a second parameter of interest is introduced as the sum of
 all design matrix entries corresponding to the positive contrast
 coefficients.  If no negative contrast coefficient is found,
 \code{pmakeDesign} sets the intercept covariates to the entries of the
 first column of the design matrix corresponding to a negative contrast coefficient.
}
\value{
  a design matrix  
}
\keyword{models}
