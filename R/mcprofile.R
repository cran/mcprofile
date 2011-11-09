mcprofile <-
function(object, CM, control=mcprofileControl(), margin=NULL, mc.cores=1) UseMethod("mcprofile")

mcprofile.glm <-
function(object, CM, control=mcprofileControl(), margin=NULL, mc.cores=1){
  if (is.null(rownames(CM))) rownames(CM) <- paste("C",1:nrow(CM), sep="")
  if (is.null(colnames(CM))) colnames(CM) <- names(coefficients(object))
  if (ncol(CM) != length(coefficients(object))) stop("Equal number of contrast and model coefficients needed!")

  df.needed <- family(object)$family == "gaussian" | length(grep("quasi", family(object)$family)) == 1 | length(grep("Negative Binomial", family(object)$family)) == 1

  control$cutoff <- qnorm(1 - control$alphamax/nrow(CM), 1)
  if (is.null(margin)){
    mmat <- NULL
  } else {
    if (nrow(rbind(margin)) == 1){
      mmat <- matrix(rep(margin[order(margin)], each=nrow(CM)), ncol=2, dimnames=list(rownames(CM), c("lower","upper")))
    } else {
      rownames(margin) <- rownames(CM)
      mmat <- margin
    }
  }

  if (mc.cores == 1){
    srdp <- list()
    optpar <- list()
    for (i in 1:nrow(CM)){
      K <- CM[i,]
      glmpro <- glm_profiling(object, K, control, margin=mmat[i,])
      srdp[[i]] <- glmpro[[1]]
      optpar[[i]] <- glmpro[[2]]
    }
  } else {
    require(parallel)
    cl <- makeCluster(mc.cores)
    glmpro <- parLapply(cl, 1:nrow(CM), function(i, CM, object, control, mmat, glm_profiling){
      K <- CM[i,]
      glm_profiling(object, K, control, margin=mmat[i,])      
    }, CM=CM, object=object, control=control, mmat=mmat, glm_profiling=mcprofile:::glm_profiling)
    stopCluster(cl)
    srdp <- lapply(glmpro, function(x) x[[1]])
    optpar <- lapply(glmpro, function(x) x[[2]])
  }
  names(srdp) <- names(optpar) <- rownames(CM)
  
  out <- list()
  out$object <- object
  out$CM <- CM
  out$srdp <- srdp
  out$optpar <- optpar
  if (df.needed) out$df <- df.residual(object) else df <- NULL
  class(out) <- "mcprofile"
  out
}

mcprofile.lm <-
function(object, CM, control=mcprofileControl(), margin=NULL, mc.cores=1){
  oc <- as.list(object$call)
  oc$family <- call("gaussian")
  oc[[1]] <- as.symbol("glm")
  object <- eval(as.call(oc))
  mcprofile.glm(object, CM=CM, control=control, margin=margin)
}


mcprofileControl <-
function(steps=8, alphamax=0.001, maxsteps=200){
  list(steps=steps, alphamax=alphamax, maxsteps=maxsteps)
}
