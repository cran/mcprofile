\name{mcpcalc}
\alias{mcpcalc}
\alias{mcpcalc.lm}
\alias{mcpcalc.glm}
\alias{mcpcalc.nls}
\alias{mcpcalc.mer}
\title{Construction of Multiple Contrast Profiles} 
\description{Calculates signed root deviance profiles given a \code{\link{glm}} or \code{\link{lm}} object.  The profiled parameters of interest are defined by providing a contrast matrix. The functions handling \code{\link{lmer}} and \code{\link{nls}} objects is only experimental; please handle with care.  }

\usage{
\method{mcpcalc}{lm}(object, CM, control=mcprofileControl(), margin=NULL, method="BFGS")
\method{mcpcalc}{glm}(object, CM, control=mcprofileControl(), margin=NULL, method="BFGS")
\method{mcpcalc}{nls}(object, CM, control=mcprofileControl(), margin=NULL, method="BFGS")
\method{mcpcalc}{mer}(object, CM, control=mcprofileControl(), margin=NULL, method="BFGS")
}

\arguments{
  \item{object}{An object of class \code{\link{glm}} or \code{\link{lm}}}
  \item{CM}{A contrast matrix with the number of columns equal the
  number of estimated parameters.  Contrast coefficients should be -1
  and 1, respectively.  Providing row names is recommendable.}
  \item{control}{A list with control arguments. See \code{\link{mcprofileControl}}.}
  \item{margin}{Numeric vector of margins, defining the profiling
  range.  If \code{NULL} no extra margin is considered.}
  \item{method}{"IRWLS", "BFGS"}
}

\details{
  According to the rows of a contrast matrix the design matrix of
a model is restructured by \code{\link{makeDesign}}. The parameter of
interest is then profiled as in function \code{\link{profile.glm}} or \code{\link{profile.nls}}.  The
signed root deviance is calculated and interpolated by interpolation splines.
}

\value{ An object inheriting from class \code{\linkS4class{mcprofile}}
  with methods \code{\link{show}}, \code{\link{plot}},
  \code{\link{confint}}, and \code{\link{test}}}

\seealso{\code{\link{profile.glm}}, \code{\link{profile.nls}}, \code{\link{glht}}, \code{\link{interpSpline}}, \code{\link{contrMat}}}


\examples{
#######################################
## cell transformation assay example ##
#######################################

str(cta)
boxplot(foci ~ conc, cta, xlab="concentration", col="lightgreen")
## change class of cta$concentration into factor
cta$concf <- factor(cta$conc, levels=unique(cta$conc))

# glm fit assuming a Poisson distribution for foci counts
# parameter estimation on the log link
# estimating marginal means by removing the intercept
fm <- glm(foci ~ concf-1, data=cta, poisson(link="log"))

### Comparing each dose to the control by Dunnett-type comparisons
# Constructing contrast matrix
CM <- contrMat(table(cta$concf), type="Dunnett")

# calculating signed root deviance profiles
(dmcp <- mcpcalc(fm, CM, method="IRWLS"))
# plot profiles
plot(dmcp)

# multiplicity adjusted p-values
(adjpv <- test(dmcp,adjust="single-step",alternative="two.sided",margin=0))
plot(adjpv, alpha=0.05, order=FALSE)
plot(dmcp, adjpv)

# simultaneous confidence intervals
(ci <- confint(dmcp, adjust="single-step"))
# exponent of confidence limits --> ratio of group means
exp(ci)

plot(exp(ci))
abline(v=1, lty=2)

plot(dmcp, ci)
}

\keyword{models}
\keyword{optimize}

