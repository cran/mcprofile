\name{mcptest-class}
\docType{class}
\alias{mcptest-class}

\title{Class "mcptest"}
\description{Objects of class "mcptest" contain slots with test
  statistics and adjusted p-values}
\section{Objects from the Class}{
  Objects are created by calls of \link[=test,mcprofile-method]{test} 
}
\section{Slots}{
  \describe{
    \item{\code{stat}:}{Numeric vector with test statistics}
    \item{\code{pvalues}:}{Numeric vector with (adjusted) p-values}
    \item{\code{margin}:}{Test margin, numeric}
    \item{\code{estimate}:}{Numeric vector of point estimates}
    \item{\code{alternative}:}{One of \code{"two.sided"}, \code{"less"}, or \code{"greater"}}
    \item{\code{adjust}:}{Multiplicity adjustment: \code{"single-step"},
    or one of the methods in \code{\link{p.adjust}}}
  }
}
\section{Methods}{
  \describe{
    \item{\link[=plot,mcptest,ANY-method]{plot}}{p-value barchart}
    \item{\link[=plot,mcprofile,mcptest-method]{plot}}{Profile plot with
      test statistics added}     
  }
}
\keyword{classes}

