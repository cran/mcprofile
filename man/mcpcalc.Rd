\name{mcpcalc}
\alias{mcpcalc}
\alias{mcpcalc.lm}
\alias{mcpcalc.glm}

\title{Construction of Multiple Contrast Profiles}

\description{Calculates signed root deviance profiles given a \code{\link{glm}} or \code{\link{lm}} object. The profiled parameters of interest are defined by providing a contrast matrix.}

\usage{
mcpcalc(object, CM, control = mcprofileControl(), margin = NULL)
\method{mcpcalc}{lm}(object, CM, control=mcprofileControl(), margin=NULL)
\method{mcpcalc}{glm}(object, CM, control=mcprofileControl(), margin=NULL)
}

\arguments{
  \item{object}{An object of class \code{\link{glm}} or \code{\link{lm}}}
  \item{CM}{A contrast matrix for the definition of parameter linear
  combinations (\code{CM \%*\% coefficients(object)}). The number of columns should be equal to the number of estimated parameters. Providing row names is recommendable.}
  \item{control}{A list with control arguments. See \code{\link{mcprofileControl}}.}
  \item{margin}{Numeric vector of margins, defining the profiling range
  with stepsize defined in the steps slot of the control argument list
  in \code{\link{mcprofileControl}}.  If \code{NULL} no extra margin is
  considered, and the stepsize is chosen automatically.}
}

\value{An object of class mcprofile. The slot \code{srdp} contains the
  profiled signed root deviance statistics.}

\details{The profiles are calculates separately for each row of the
  contrast matrix. With single non-zero entries in a contrast matrix,
  the \code{\link{glm.fit}} function is used for optimization. With two
  non-zero entries in a contrast, one-dimensional optimization is
  performed by \code{\link{optimize}} for the second non-zero contrast parameter conditional on
  the parameter linear combination and all other nuisance
  parameters. For all higher dimensional optimization problems, the BFGS
  method in \code{\link{optim}} is used.
}

\author{Daniel Gerhard}

\seealso{\code{\link[MASS]{profile.glm}}, \code{\link[multcomp]{glht}},
  \code{\link[multcomp]{contrMat}}, \code{\link{confint.mcprofile}}, \code{\link{summary.mcprofile}}}

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
(dmcp <- mcpcalc(fm, CM))
# plot profiles
plot(dmcp)
# confidence intervals
(ci <- confint(dmcp))
plot(ci)
}

\keyword{ models }
\keyword{ optimize }
