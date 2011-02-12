\name{mcprofileRatio-class}
\docType{class}
\alias{mcprofileRatio-class}
\alias{initialize, mcprofileRatio-method}

\title{Class "mcprofileRatio"}
\description{Objects of class "mcprofileRatio" contain signed root deviance
  profiles and some corresponding interpolation splines}
\section{Objects from the Class mcprofile}{
  Objects are created by calls of \code{\link{mcpcalcRatio}}
}
\section{Slots}{
  \describe{
    \item{\code{CMn}:}{Numerator contrast matrix}
    \item{\code{CMd}:}{Denominator contrast matrix}
    \item{\code{control}:}{control arguments of class \code{"list"} with
      slots steps, alphamax, and fixed.range}
    \item{\code{method}:}{One of \code{"IRWLS"}, \code{"BFGS"}, or
      \code{"constrOptim"} depending on the optimization method which has been used}    
    \item{\code{estimate}:}{Numeric vector of parameter estimates} 
    \item{\code{model}:}{Object of class "list" with information aout
      the fitted model, see \code{\link{glmobj}}}
    \item{\code{SRDP}:}{Signed root deviance profile object of class
      \code{"list"}, for each row of the contrast matrix a slot is
      available, containing a dataframe with the signed root
      deviance estimates z, the parameter of interest b, and the
      conditional variance estimates}
    \item{\code{fsplines}:}{A list of interpolation spline functions predicting z given b}
  }
}
\section{Methods}{
  \describe{
    \item{\link[=plot,mcprofileRatio,ANY-method]{plot}}{A lattice profile plot}
    \item{\link{test}}{p-values for general linear hypotheses tests}
    \item{\link[=confint,mcprofileRatio-method]{confint}}{Simultaneous confidence intervals}
  }
}
\seealso{\code{\link{mcpcalcRatio}}}
\keyword{classes}


