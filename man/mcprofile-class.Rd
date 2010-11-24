\name{mcprofile-class}
\docType{class}
\alias{mcprofile-class}
\alias{initialize, mcprofile-method}

\title{Class "mcprofile"}
\description{Objects of class "mcprofile" contain signed root deviance
  profiles and some corresponding interpolation splines}
\section{Objects from the Class mcprofile}{
  Objects are created by calls of \code{\link{mcpcalc}}
}
\section{Slots}{
  \describe{
    \item{\code{CM}:}{Contrast matrix}
    \item{\code{control}:}{control arguments of class \code{"list"} with
      slots steps, alphamax, and fixed.range}
    \item{\code{method}:}{One of \code{"IRWLS"}, \code{"BFGS"}, or
      \code{"constrOptim"} depending on the optimization method which has been used}    
    \item{\code{estimate}:}{Numeric vector of parameter estimates} 
    \item{\code{vest}:}{Numeric vector of scaled variance estimates}
    \item{\code{vestunsc}:}{Numeric vector of unscaled variance estimates}
    \item{\code{dvest}:}{Numeric vector of unscaled variance estimates
      of the original model}    
    \item{\code{model}:}{Object of class "list" with information aout
      the fitted model, see \code{\link{glmobj}}}
    \item{\code{SRDP}:}{Signed root deviance profile object of class
      \code{"list"}, for each row of the contrast matrix a slot is
      available, containing a dataframe with the signed root
      deviance estimates z, the parameter of interest b, and the
      conditional variance estimates}
    \item{\code{fsplines}:}{A list of interpolation spline functions predicting z given b}
    \item{\code{bsplines}:}{A list of interpolation spline functions predicting b given z}
  }
}
\section{Methods}{
  \describe{
    \item{\link[=plot,mcprofile,ANY-method]{plot}}{A lattice profile plot}
    \item{\link{wald}}{Approximate profiles by a quadratic function}
    \item{\link{hoa}}{Higher order approximations}
    \item{\link{test}}{p-values for general linear hypotheses tests}
    \item{\link[=confint,mcprofile-method]{confint}}{Simultaneous confidence intervals}
  }
}
\seealso{\code{\link{mcpcalc}}}
\keyword{classes}


