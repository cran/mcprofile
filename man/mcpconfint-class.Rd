\name{mcpconfint-class}
\docType{class}
\alias{mcpconfint-class}

\title{Class "mcpconfint"}
\description{Objects of class "mcpconfint" contain slots with estimated
  confidence limits}
\section{Objects from the Class}{
  Objects are created by calls of \link[=confint,mcprofile-method]{confint}
}
\section{Slots}{
  \describe{
    \item{\code{confint}:}{A dataframe with lower and/or upper confidence limits}
    \item{\code{quantile}:}{The signed root deviance cutoff value}
    \item{\code{estimate}:}{A numeric vector with point estimates}
    \item{\code{alternative}:}{One of \code{"two.sided"}, \code{"less"}, or \code{"greater"}} 
    \item{\code{adjust}:}{Multiplicity adjustment: \code{"none"}, \code{"bonferroni"}, or \code{"single-step"}}
  }
}
\section{Methods}{
  \describe{
    \item{\link[=plot,mcpconfint,ANY-method]{plot}}{Confidence interval plot}
    \item{\link[=plot,mcprofile,mcpconfint-method]{plot}}{Profile plot with confidence limits added}
  }
}
\keyword{classes}

