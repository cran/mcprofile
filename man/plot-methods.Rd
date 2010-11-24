\name{plot-methods}
\docType{methods}
\alias{plot,mcpconfint,ANY-method}
\alias{plot,mcprofile,ANY-method}
\alias{plot,mcprofile,mcprofile-method}
\alias{plot,mcprofile,mcpconfint-method}
\alias{plot,mcprofile,mcptest-method}
\alias{plot,mcptest,ANY-method}
\alias{plot,mcprofileRatio,ANY-method}
\alias{plot,mcprofileRatio,mcpconfint-method}
\alias{plot,mcprofileRatio,mcptest-method}
\title{Methods for Function plot in Package `mcprofile'}
\description{
 Display multiple contrast profiles, confidence intervals, and p-values.
}

\usage{
\S4method{plot}{mcprofile,ANY}(x, y, squared=FALSE, trans=identity, ...)
\S4method{plot}{mcpconfint,ANY}(x, y, margin=NULL, pch=19, trans=identity, ...)
\S4method{plot}{mcptest,ANY}(x, y, trans=identity, alpha=NULL, order=TRUE, ...)
\S4method{plot}{mcprofile,mcpconfint}(x, y, squared=FALSE, trans=identity, ...)
\S4method{plot}{mcprofile,mcptest}(x, y, squared=FALSE, trans=identity, ...)
\S4method{plot}{mcprofile,mcprofile}(x, y, squared=FALSE,
trans=identity, ...)
\S4method{plot}{mcprofileRatio,ANY}(x, y, squared=FALSE, trans=identity,
...)
\S4method{plot}{mcprofileRatio,mcpconfint}(x, y, squared=FALSE, trans=identity, ...)
\S4method{plot}{mcprofileRatio,mcptest}(x, y, squared=FALSE, trans=identity, ...)
}

\arguments{
  \item{x}{An object of class \code{\linkS4class{mcprofile}}, or \code{\linkS4class{mcprofileRatio}}}
  \item{y}{Missing, or an object of class \code{\linkS4class{mcpconfint}}, or \code{\linkS4class{mcptest}}}
  \item{squared}{Take the square of the signed root deviance}
  \item{trans}{Function to transform the x-axis, default \code{\link{identity}}}
  \item{margin}{Include values to the x-range of the confidence intervals plot}
  \item{pch}{Point character for point estimates in confidence interval plot}
  \item{alpha}{Add line at value of alpha in p-value plot}
  \item{order}{Order p-value bars by size}
  \item{\dots}{optional arguments to \code{\link{xyplot}}, or \code{\link{plot}}}
}

\section{Methods}{
  \describe{
    \item{\code{plot(x="mcprofile")}}{Profile plot}
    \item{\code{plot(x="mcprofileRatio")}}{Profile plot for ratios of parameters}
    \item{\code{plot(x="mcpconfint")}}{Confidence interval plot}
    \item{\code{plot(x="mcptest")}}{p-value plot}
    \item{\code{plot(x="mcprofile", y="mcpconfint")}}{Add confidence
      limits to profile plots}
    \item{\code{plot(x="mcprofileRatio", y="mcpconfint")}}{Add confidence
      limits to a ratio profile plot}
    \item{\code{plot(x="mcprofile", y="mcptest")}}{Add test
      statistics to profile plots}
    \item{\code{plot(x="mcprofileRatio", y="mcptest")}}{Add test
      statistics to a ratio profile plot}
    \item{\code{plot(x="mcprofile", y="mcprofile")}}{Compare two profile
  objects}
  }
}
\details{
  Function \code{\link{xyplot}} is used to display the profiles. 
}
\keyword{methods}
\keyword{hplot}
