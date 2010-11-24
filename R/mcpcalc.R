mcpcalc <- function(object, CM, control=mcprofileControl(), margin=NULL, method="BFGS") UseMethod("mcpcalc")

mcpcalc.lm <- function(object, CM, control=mcprofileControl(), margin=NULL, method="BFGS"){
  object <- glm(formula(object), data=object$model, family=gaussian(link="identity"))
  if (method == "IRWLS") return(mcpcalcIRWLS(object, CM, control, margin))
  if (method == "BFGS") return(mcpcalcBFGS(object, CM, control, margin))
}

mcpcalc.glm <- function(object, CM, control=mcprofileControl(), margin=NULL, method="BFGS"){
  if (method == "IRWLS") return(mcpcalcIRWLS(object, CM, control, margin))
  if (method == "BFGS") return(mcpcalcBFGS(object, CM, control, margin))
}

mcpcalc.nls <- function(object, CM, control=mcprofileControl(), margin=NULL, method="BFGS"){
  if (method == "IRWLS") stop("IRWLS optimisation not available.")
  if (method == "BFGS") return(mcpcalcnls(object, CM, control, margin))      
}
#############################

mcpcalcIRWLS <- function(object, CM, control=mcprofileControl(), margin=NULL){
  if (is.null(rownames(CM))) rownames(CM) <- paste("C",1:nrow(CM), sep="")
  if (is.null(colnames(CM))) colnames(CM) <- names(coefficients(object))
  if (ncol(CM) != length(coefficients(object))) stop("Equal number of contrast and model coefficients needed!")
  if (control$fixed.range & length(margin) < 2) stop("2 margins have to be supplied to fix the step sizes between them (fixed.range = TRUE)!")
  cutoff <- qnorm(1 - control$alphamax/nrow(CM), 1)
  delta <- cutoff/control$steps
  mod <- glmobj(object)
  off <- mod$offset   
  SRDP <- list()
  estimate <- numeric(length=nrow(CM))
  vest <- numeric(length=nrow(CM))
  vestunsc <- numeric(length=nrow(CM))
  dvest <- numeric(length=nrow(CM))
  for (i in 1:nrow(CM)){
    K <- CM[i,]
    if (any(K < 0) & any(K > 0)){
      Xi <- makeDesign(K, mod$X)
      startc <- c(mod$coefficients[1], K %*% mod$coefficients, mod$coefficients[K == 0])
    }
    if (all(K >= 0)){
      Xi <- pmakeDesign(K, mod$X)
      startc <- c(mod$coefficients[K == 0][1], K %*% mod$coefficients, mod$coefficients[K == 0][-1])
    }
    pfm <- ofm <- condOpt(mod, Xi, off, startc)
    Xio <- Xi[,-2,drop=FALSE]
    Xib <- Xi[,2,drop=FALSE]
    estimate[i] <- pp1 <- ofm$coefficients[2]
    vest[i] <- ofm$vcov[2,2]
    vestunsc[i] <- ofm$cov.unscaled[2,2]
    dvest[i] <- det(ofm$cov.unscaled)
    if (is.null(margin)) minmargin <- pp1 else  minmargin <- min(margin)
    if (is.null(margin)) maxmargin <- pp1 else  maxmargin <- max(margin)
    if (pp1 != 0) pp2 <- pp1+abs(pp1)*0.01 else pp2 <- 0.01
    if (control$fixed.range){
      stepsize <- diff(c(pp1, max(margin))) / control$steps
      pp2 <- pp1 + stepsize
    }
    zp1 <- 0
    z <- 0
    b <- pp1
    dv <- det(ofm$cov.unscaled)
    startc <- pfm$coefficients[-2]
    mst <- 1
    while((zp1 < cutoff | pp1 < maxmargin) & pfm$converged & mst <= control$maxsteps){
      offp <- off + (Xib %*% pp2)[,1]
      pfm <- try(condOpt(pfm, Xio, offp, startc), silent=TRUE)
      if (class(pfm) == "try-error") break
      if (pfm$converged){
        zp <- srdp(pfm, ofm)
        if (!control$fixed.range) stepsize <- sign(pp2-pp1)*(delta*((pp2-pp1)/(zp-zp1)))
        pp <- pp2 + stepsize
        pp1 <- pp2
        pp2 <- pp
        zp1 <- zp
        z <- c(z, zp1)
        b <- c(b, pp1)
        dv <- c(dv, det(pfm$cov.unscaled))
        startc <- pfm$coefficients
        if (control$fixed.range & length(b) == control$steps+1) break
        mst <- mst + 1
      }
    }
    pfm <- ofm
    pp1 <- ofm$coefficients[2]
    if (pp1 != 0) pp2 <- pp1-abs(pp1*0.01) else pp2 <- -0.01
    if (control$fixed.range){
      stepsize <- diff(c(min(margin), pp1)) / control$steps
      pp2 <- pp1 - stepsize
    }    
    zp1 <- 0
    startc <- pfm$coefficients[-2]
    mst <- 1
    while((zp1 < cutoff | pp1 > minmargin) & pfm$converged & mst <= control$maxsteps){
      offp <- off + (Xib %*% pp2)[,1]
      pfm <- try(condOpt(pfm, Xio, offp, startc), silent=TRUE)
      if (class(pfm) == "try-error") break
      if (pfm$converged){
        zp <- srdp(pfm, ofm)
        if (!control$fixed.range) stepsize <- sign(pp2-pp1)*(delta*((pp2-pp1)/(zp-zp1)))
        pp <- pp2 - stepsize
        pp1 <- pp2
        pp2 <- pp
        zp1 <- zp
        z <- c(z, -zp1)
        b <- c(b, pp1)
        dv <- c(dv, det(pfm$cov.unscaled))
        startc <- pfm$coefficients
        if (control$fixed.range & length(b) == 2*control$steps+1) break
        mst <- mst + 1
      }
    }
    zd <- data.frame(z, b, dv)
    mat <- zd[order(zd[,2]),]
    mat <- unique(na.omit(mat[!is.nan(mat[,1]),]))
    SRDP[[i]] <- mat
  }
  names(estimate) <- names(SRDP) <- rownames(CM)
  ########
  fsplist <- lapply(SRDP, function(z){
    try(interpSpline(z[,2], z[,1]), silent=TRUE)
  })
  bsplist <- lapply(SRDP, function(z){
    try(interpSpline(z[,1], z[,2]), silent=TRUE)
  })
  new(Class="mcprofile", CM=CM, estimate=estimate, vest=vest, vestunsc=vestunsc, dvest=dvest, model=mod, SRDP=SRDP, fsplines=fsplist, bsplines=bsplist, control=control, method="IRWLS")
}


#####################################################
#####################################################
#####################################################


mcpcalcBFGS <- function(object, CM, control=mcprofileControl(), margin=NULL){
  if (is.null(rownames(CM))) rownames(CM) <- paste("C",1:nrow(CM), sep="")
  if (is.null(colnames(CM))) colnames(CM) <- names(coefficients(object))
  if (ncol(CM) != length(coefficients(object))) stop("Equal number of contrast and model coefficients needed!")
  if (control$fixed.range & length(margin) < 2) stop("2 margins have to be supplied to fix the step sizes between them (fixed.range = TRUE)!")
  kbb <- function(x, k, p, w){
    K <- -1*c(-p, k[-w])/k[w]
    hv <- c(1,x)
    b1 <- (K %*% hv)[1,1]
    b <- numeric(length=length(k))
    b[(1:length(k))[-w]] <- x
    b[w] <- b1
    b
  } 
  cutoff <- qnorm(1 - control$alphamax/nrow(CM), 1)
  delta <- cutoff/control$steps
  mod <- glmobj(object)
  off <- mod$offset
  X <- mod$X
  y <- mod$response
  dev0 <- mod$deviance
  disp <- mod$dispersion
  est <- coefficients(object)
  estimate <- as.vector(CM %*% coefficients(object))
  vest <- diag(CM %*% mod$vcov %*% t(CM))
  vestunsc <- diag(CM %*% mod$cov.unscaled %*% t(CM))
  dvest <- numeric(length=length(estimate))
  .denvi <- new.env()
  .start <- new.env()
  SRDP <- list()
  for (i in 1:nrow(CM)){
    K <- CM[i,]
    pos <- abs(K) > 0
    neg <- abs(K) == 0
    b1 <- b <- as.vector(K %*% est)
    if (is.null(margin)) minmargin <- b1 else  minmargin <- min(margin)
    if (is.null(margin)) maxmargin <- b1 else  maxmargin <- max(margin)
    if (b1 != 0) b2 <- b1+abs(b1)*0.01 else b2 <- 0.01
    if (control$fixed.range){
      stepsize <- diff(c(b1, max(margin))) / control$steps
      b2 <- b1 + stepsize
    }
    z <- z1 <- z2 <- 0
    w <- which(K != 0)[1]
    xi <- X[,neg, drop=FALSE]
    bp <- est[pos]
    dv <- det(mod$cov.unscaled)
    bstart <- est[-w]
    assign("startc", est[neg], envir = .start)
    mst <- 1
    while ((z2 < cutoff | b1 < maxmargin) & mst <= control$maxsteps){
      optd <- function(x){
        Ki <- -1*c(-b2, K[-w])/K[w]
        hv <- c(1,x)
        bo <- (Ki %*% hv)[1,1]
        bopt <- numeric(length=length(K))
        bopt[(1:length(K))[-w]] <- x
        bopt[w] <- bo
        offs <- X[,!neg,drop=FALSE] %*% bopt[!neg] + off
        startc <- get("startc", envir=.start)
        cfm <- glm.fit(xi,y,weights=mod$weights,etastart=mod$fitted,offset=offs,family=mod$family,control=mod$control, start=startc)
        class(cfm) <- "glm"
        assign(".DV", summary(cfm)$cov.unscaled, envir = .denvi)
        assign("startc", cfm$coefficients, envir = .start)
        cfm$deviance
      }
      conopt <- try(optim(bstart, optd, method="BFGS"), silent=TRUE)
      if (class(conopt) == "try-error") break
      opar <- kbb(conopt$par, K, b2, w)
      bstart <- opar[-w]
      dev1 <- conopt$value
      z2 <- sqrt(abs(dev1 - dev0)/disp)
      z <- c(z, z2)
      b <- c(b, b2)
      dv <- c(dv, det(get(".DV", envir=.denvi)))
      if (!control$fixed.range) stepsize <- sign(b2-b1)*(delta*((b2-b1)/(z2-z1)))
      b1 <- b2
      b2 <- b2 + stepsize
      z1 <- z2
      if (control$fixed.range & length(b) == 2*control$steps+1) break
      mst <- mst + 1
    }
    b1 <- as.vector(K %*% est)
    if (b1 != 0) b2 <- b1-abs(b1)*0.01 else b2 <- -0.01
    if (control$fixed.range){
      stepsize <- diff(c(min(margin), b1)) / control$steps
      b2 <- b1 - stepsize
    }
    z1 <- z2 <- 0
    bstart <- est[-w]
    assign("startc", est[neg], envir = .start)
    mst <- 1
    while ((z2 < cutoff | b1 > minmargin) & mst <= control$maxsteps){
      optd <- function(x){
        Ki <- -1*c(-b2, K[-w])/K[w]
        hv <- c(1,x)
        bo <- (Ki %*% hv)[1,1]
        bopt <- numeric(length=length(K))
        bopt[(1:length(K))[-w]] <- x
        bopt[w] <- bo
        offs <- X[,!neg,drop=FALSE] %*% bopt[!neg] + off
        startc <- get("startc", envir=.start)
        cfm <- glm.fit(xi,y,weights=mod$weights,etastart=mod$fitted,offset=offs,family=mod$family,control=mod$control, start=startc)
        class(cfm) <- "glm"
        assign(".DV", summary(cfm)$cov.unscaled, envir = .denvi)
        assign("startc", cfm$coefficients, envir = .start)
        cfm$deviance
      }
      conopt <- try(optim(bstart, optd, method="BFGS"), silent=TRUE)
      if (class(conopt) == "try-error") break
      opar <- kbb(conopt$par, K, b2, w)
      bstart <- opar[-w]      
      dev1 <- conopt$value
      z2 <- sqrt(abs(dev1 - dev0)/disp)
      z <- c(z, -z2)
      b <- c(b, b2)
      dv <- c(dv, det(get(".DV", envir=.denvi)))
      if (!control$fixed.range) stepsize <- sign(b2-b1)*(delta*((b2-b1)/(z2-z1)))
      b1 <- b2
      b2 <- b2 - stepsize
      z1 <- z2
      if (control$fixed.range & length(b) == 2*control$steps+1) break
      mst <- mst + 1
    }
    zd <- data.frame(z, b, dv)
    mat <- zd[order(zd[,2]),]
    mat <- unique(na.omit(mat[!is.nan(mat[,1]),]))
    SRDP[[i]] <- mat
    nuiKv <- diag(neg)
    Kv <- rbind(K, nuiKv[apply(nuiKv,1,sum) != 0,])
    dvest[i] <- det(Kv %*% mod$cov.unscaled %*% t(Kv))
  }
  names(SRDP) <- rownames(CM)
  ########
  fsplist <- lapply(SRDP, function(z){
    try(interpSpline(z[,2], z[,1]), silent=TRUE)
  })
  bsplist <- lapply(SRDP, function(z){
    try(interpSpline(z[,1], z[,2]), silent=TRUE)
  })
  new(Class="mcprofile", CM=CM, estimate=estimate, vest=vest, vestunsc=vestunsc, dvest=dvest, model=mod, SRDP=SRDP, fsplines=fsplist, bsplines=bsplist, control=control, method="BFGS")
}





##########################################################################

mcpcalcnls <- function(object, CM, control=mcprofileControl(), margin=NULL){
  if (is.null(rownames(CM))) rownames(CM) <- paste("C",1:nrow(CM), sep="")
  if (is.null(colnames(CM))) colnames(CM) <- names(coefficients(object))
  if (ncol(CM) != length(coefficients(object))) stop("Equal number of contrast and model coefficients needed!")
  if (control$fixed.range & length(margin) < 2) stop("2 margins have to be supplied to fix the step sizes between them (fixed.range = TRUE)!")  
  kbb <- function(x, k, p, w){
    K <- -1*c(-p, k[-w])/k[w]
    hv <- c(1,x)
    b1 <- (K %*% hv)[1,1]
    b <- numeric(length=length(k))
    b[(1:length(k))[-w]] <- x
    b[w] <- b1
    b
  } 
  cutoff <- qnorm(1 - control$alphamax/nrow(CM), 1)
  delta <- cutoff/control$steps
  prof <- stats:::profiler.nls(object)
  pars <- est <- prof$getFittedPars()
  estimate <- as.vector(CM %*% coefficients(object))
  vest <- diag(CM %*% vcov(object) %*% t(CM))
  dvest <- numeric(length=length(estimate))
  SRDP <- list()
  for (i in 1:nrow(CM)){
    k <- CM[i,]
    pos <- abs(k) > 0
    neg <- abs(k) == 0
    b1 <- b <- as.vector(k %*% pars)
    if (is.null(margin)) minmargin <- b1 else  minmargin <- min(margin)
    if (is.null(margin)) maxmargin <- b1 else  maxmargin <- max(margin)    
    K <- k[!neg]
    if (b1 != 0) b2 <- b1+abs(b1)*0.01 else b2 <- 0.01
    if (control$fixed.range){
      stepsize <- diff(c(b1, max(margin))) / control$steps
      b2 <- b1 + stepsize
    }
    z <- z1 <- z2 <- 0
    w <- which(K != 0)[1]
    bp <- est[pos]
    bstart <- est[!neg][-w]
    par <- which(!neg)
    pars <- prof$getFittedPars()
    prof$setDefault(varying = par)
    mst <- 1
    while (z2 < cutoff & mst <= control$maxsteps){
      optfunc <- function(x){
        Ki <- -1*c(-b2, K[-w])/K[w]
        hv <- c(1,x)
        bo <- (Ki %*% hv)[1,1]
        bopt <- numeric(length=length(K))
        bopt[(1:length(K))[-w]] <- x
        bopt[w] <- bo
        pars[par] <- bopt
        prof$setDefault(params = pars)
        prof$setDefault(varying = par)
        ans <- prof$getProfile()
        fst <- ans$fstat
        prof$setDefault()
        fst
      }
      conopt <- try(optim(bstart, optfunc, method="BFGS"), silent=TRUE)
      if (class(conopt) == "try-error") break
      opar <- kbb(conopt$par, K, b2, w)
      bstart <- opar[-w]
      z2 <- sqrt(conopt$value)
      z <- c(z, z2)
      b <- c(b, b2)
      if (!control$fixed.range) stepsize <- sign(b2-b1)*(delta*((b2-b1)/(z2-z1)))
      b1 <- b2
      b2 <- b2 + stepsize
      z1 <- z2
      if (control$fixed.range & length(b) == 2*control$steps+1) break
      mst <- mst + 1
    }
    b1 <- as.vector(k %*% pars)
    if (b1 != 0) b2 <- b1-abs(b1)*0.01 else b2 <- -0.01
    if (control$fixed.range){
      stepsize <- diff(c(min(margin), b1)) / control$steps
      b2 <- b1 - stepsize
    }
    z1 <- z2 <- 0
    bstart <- est[!neg][-w]
    mst <- 1
    while (z2 < cutoff & mst <= control$maxsteps){
      optfunc <- function(x){
        Ki <- -1*c(-b2, K[-w])/K[w]
        hv <- c(1,x)
        bo <- (Ki %*% hv)[1,1]
        bopt <- numeric(length=length(K))
        bopt[(1:length(K))[-w]] <- x
        bopt[w] <- bo
        pars[par] <- bopt
        prof$setDefault(params = pars)
        prof$setDefault(varying = par)
        ans <- prof$getProfile()
        fst <- ans$fstat
        prof$setDefault()
        fst
      }
      conopt <- try(optim(bstart, optfunc, method="BFGS"), silent=TRUE)
      if (class(conopt) == "try-error") break
      opar <- kbb(conopt$par, K, b2, w)
      bstart <- opar[-w]
      z2 <- sqrt(conopt$value)
      z <- c(z, -z2)
      b <- c(b, b2)
      if (!control$fixed.range) stepsize <- sign(b2-b1)*(delta*((b2-b1)/(z2-z1)))
      b1 <- b2
      b2 <- b2 - stepsize
      z1 <- z2
      if (control$fixed.range & length(b) == 2*control$steps+1) break
      mst <- mst + 1
    }
    zd <- data.frame(z, b)
    mat <- zd[order(zd[,2]),]
    mat <- unique(na.omit(mat[!is.nan(mat[,1]),]))
    SRDP[[i]] <- mat
    nuiKv <- diag(neg)
    Kv <- rbind(k, nuiKv[apply(nuiKv,1,sum) != 0,])
    dvest[i] <- det(Kv %*% vcov(object) %*% t(Kv))
  }
  names(SRDP) <- rownames(CM)
  ########
  fsplist <- lapply(SRDP, function(z){
    try(interpSpline(z[,2], z[,1]), silent=TRUE)
  })
  bsplist <- lapply(SRDP, function(z){
    try(interpSpline(z[,1], z[,2]), silent=TRUE)
  })
  new(Class="mcprofile", CM=CM, estimate=estimate, vest=vest, vestunsc=vest, dvest=dvest, model=list(df=df.residual(object)), SRDP=SRDP, fsplines=fsplist, bsplines=bsplist, control=control, method="BFGS")
}
  


