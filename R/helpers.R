makeDesign <- function(k, X){
  pos <- k > 0
  neg <- k < 0
  zero <- k == 0
  Xx <- X[,!zero,drop=FALSE]
  Xp <- X[,pos,drop=FALSE]
  Xint <- apply(Xx, 1, sum)
  Xpar <- apply(Xp, 1, sum)
  Xnui <- X[,zero,drop=FALSE]
  cbind(Xint, Xpar, Xnui)
}

pmakeDesign <- function(k, X){
  pos <- k > 0
  zero <- k == 0
  wz <- which(zero)[1]
  zero[wz] <- FALSE
  Xx <- X[,wz,drop=FALSE]
  Xp <- X[,pos,drop=FALSE]
  Xint <- apply(Xx, 1, sum)
  Xpar <- apply(Xp, 1, sum)
  Xnui <- X[,zero,drop=FALSE]
  cbind(Xint, Xpar, Xnui)
}

glmobj <- function(object){
  off <- if (is.null(model.offset(object))) rep(0, length(object$y)) else model.offset(object)
	dfneeded <- family(object)$family == "gaussian" | length(grep("quasi", family(object)$family)) == 1 | length(grep("Negative Binomial", family(object)$family)) == 1
  return(list(coefficients=coefficients(object),
       vcov=vcov(object),
       cov.unscaled=summary(object)$cov.unscaled,
       X=model.matrix(object),
       deviance=deviance(object),
       dispersion=summary(object)$dispersion,
       response=object$y,
       weights=weights(object),
       family=family(object),
       fitted=fitted(object),
       control=object$control,
       offset=off,
       df=if (dfneeded) df.residual(object) else NA,
       converged=object$converged))
}



## conditional optimization of glms
condOpt <- function(object, X, offset, start=NULL){
  object$offset <- offset
  object$X <- X
  fm <- glm.fit(object$X, object$response, weights=object$weights,etastart=object$fitted,offset=object$offset,family=object$family,control=glm.control(maxit = 2000), start=start)
  class(fm) <- c("glm")
  object$coefficients <- coefficients(fm)
  object$deviance <- deviance(fm)
  object$fitted <- fitted(fm)
  object$vcov <- vcov(fm)
  object$cov.unscaled <- summary(fm)$cov.unscaled
  object$converged <- fm$converged
  return(object)
}



### deviance calculation
srdp <- function(cobject, oobject){
  devdiff <- abs(cobject$deviance - oobject$deviance)
  stat <- sqrt(devdiff/oobject$dispersion)
  return(stat)
}


### expit
expit <- function(x) exp(x)/(1 + exp(x))


mcprofileControl <- function(steps=8, alphamax=0.001, fixed.range=FALSE, maxsteps=200){
  list(steps=steps, alphamax=alphamax, fixed.range=fixed.range, maxsteps=maxsteps)
}
