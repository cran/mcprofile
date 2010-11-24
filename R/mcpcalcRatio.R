mcpcalcRatio <- function(object, CMn, CMd, control=mcprofileControl(), margin=NULL, method="BFGS") UseMethod("mcpcalcRatio")

mcpcalcRatio.lm <- function(object, CMn, CMd, control=mcprofileControl(), margin=NULL, method="BFGS"){
  object <- glm(formula(object), data=object$model, family=gaussian(link="identity"))
  if (method == "IRWLS") return("Not yet implemented ...")
  if (method == "BFGS") return(mcpcalcRatioBFGS(object, CMn, CMd, control, margin))
}

mcpcalcRatio.glm <- function(object, CMn, CMd, control=mcprofileControl(), margin=NULL, method="BFGS"){
  if (method == "IRWLS") return("Not yet implemented ...")
  if (method == "BFGS") return(mcpcalcRatioBFGS(object, CMn, CMd, control, margin))
}

mcpcalcRatio.nls <- function(object, CMn, CMd, control=mcprofileControl(), margin=NULL, method="BFGS"){
  if (method == "IRWLS") stop("IRWLS optimisation not available.")
  if (method == "BFGS") return(mcpcalcnlsRatio(object, CMn, CMd, control, margin))      
}
#############################

#####################################################
#####################################################


mcpcalcRatioBFGS <- function(object, CMn, CMd, control=mcprofileControl(), margin=NULL){
  if (is.null(rownames(CMn))) rownames(CMn) <- paste("C",1:nrow(CMn), sep="")
  if (is.null(colnames(CMn))) colnames(CMn) <- names(coefficients(object))
  if (ncol(CMn) != length(coefficients(object)) | ncol(CMd) != length(coefficients(object))) stop("Equal number of contrast and model coefficients needed!")
  if (control$fixed.range & length(margin) < 2) stop("2 margins have to be supplied to fix the step sizes between them (fixed.range = TRUE)!")
  kbbr <- function(x, Kn,Kd, b2, w){
    bw <- (b2*Kd[-w] %*% x - Kn[-w] %*% x) / (Kn[w] - b2*Kd[w])
    b <- numeric(length=length(Kn))
    b[(1:length(Kn))[-w]] <- x
    b[w] <- bw
    b
  }
  cutoff <- qnorm(1 - control$alphamax/nrow(CMn), 1)
  delta <- cutoff/control$steps
  mod <- glmobj(object)
  off <- mod$offset
  X <- mod$X
  y <- mod$response
  dev0 <- mod$deviance
  disp <- mod$dispersion
  est <- coefficients(object)
  estimate <- as.vector(CMn %*% est / CMd %*% est)
  .start <- new.env()
  SRDP <- list()
  for (i in 1:nrow(CMn)){
    Kn <- CMn[i,]
    Kd <- CMd[i,]
    neg <- Kn == 0 & Kd == 0
    w <- which(Kn != 0)[1]
    bstart <- est[-w]
    b <- b1 <- (Kn %*% est / Kd %*% est)[1,1]    
    if (is.null(margin)) minmargin <- b1 else  minmargin <- min(margin)
    if (is.null(margin)) maxmargin <- b1 else  maxmargin <- max(margin)
    if (b1 != 0) b2 <- b1+abs(b1)*0.01 else b2 <- 0.01
    if (control$fixed.range){
      stepsize <- diff(c(b1, max(margin))) / control$steps
      b2 <- b1 + stepsize
    }
    z <- z1 <- z2 <- 0
    xi <- X[,neg, drop=FALSE]
    assign("startc", est[neg], envir = .start)
    mst <- 1
    while (z2 < cutoff & mst <= control$maxsteps){
      optd <- function(x){
        bw <- (b2*Kd[-w] %*% x - Kn[-w] %*% x) / (Kn[w] - b2*Kd[w])
        bopt <- numeric(length=length(Kn))
        bopt[(1:length(Kn))[-w]] <- x
        bopt[w] <- bw
        offs <- X[,!neg, drop=FALSE] %*% bopt[!neg] + off
        startc <- get("startc", envir=.start)
        cfm <- glm.fit(xi,y,weights=mod$weights,etastart=mod$fitted,offset=offs,family=mod$family,control=mod$control)
        assign("startc", cfm$coefficients, envir = .start)
        cfm$deviance
      }
      conopt <- try(optim(bstart, optd, method="BFGS"), silent=TRUE)
      if (class(conopt) == "try-error") break
      opar <- kbbr(conopt$par, Kn,Kd, b2, w)
      bstart <- opar[-w]
      dev1 <- conopt$value
      z2 <- sqrt(abs(dev1 - dev0)/disp)
      z <- c(z, z2)
      b <- c(b, b2)
      if (!control$fixed.range) stepsize <- sign(b2-b1)*(delta*((b2-b1)/(z2-z1)))
      b1 <- b2
      b2 <- b2 + stepsize
      z1 <- z2
      if (control$fixed.range & length(b) == 2*control$steps+1) break
      mst <- mst + 1
    }
    b1 <- (Kn %*% est / Kd %*% est)[1,1]
    if (b1 != 0) b2 <- b1-abs(b1)*0.01 else b2 <- -0.01
    if (control$fixed.range){
      stepsize <- diff(c(min(margin), b1)) / control$steps
      b2 <- b1 - stepsize
    }
    z1 <- z2 <- 0
    bstart <- est[-w]
    assign("startc", est[neg], envir = .start)
    mst <- 1
    while (z2 < cutoff & mst <= control$maxsteps){
      optd <- function(x){
        bw <- (b2*Kd[-w] %*% x - Kn[-w] %*% x) / (Kn[w] - b2*Kd[w])
        bopt <- numeric(length=length(Kn))
        bopt[(1:length(Kn))[-w]] <- x
        bopt[w] <- bw
        offs <- X[,!neg, drop=FALSE] %*% bopt[!neg] + off
        startc <- get("startc", envir=.start)
        cfm <- glm.fit(xi,y,weights=mod$weights,etastart=mod$fitted,offset=offs,family=mod$family,control=mod$control)
        assign("startc", cfm$coefficients, envir = .start)
        cfm$deviance
      }
      conopt <- try(optim(bstart, optd, method="BFGS"), silent=TRUE)
      if (class(conopt) == "try-error") break
      opar <- kbbr(conopt$par, Kn,Kd, b2, w)
      bstart <- opar[-w]      
      dev1 <- conopt$value
      z2 <- sqrt(abs(dev1 - dev0)/disp)
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
  }
  names(estimate) <- names(SRDP) <- rownames(CMn)
  ########
  fsplist <- lapply(SRDP, function(z){
    try(interpSpline(z[,2], z[,1]), silent=TRUE)
  })
  bsplist <- lapply(SRDP, function(z){
    try(interpSpline(z[,1], z[,2]), silent=TRUE)
  })
  new(Class="mcprofileRatio", CMn=CMn, CMd=CMd, estimate=estimate, model=mod, SRDP=SRDP, fsplines=fsplist, bsplines=bsplist, control=control, method="BFGS")
}



mcpcalcnlsRatio <- function(object, CMn, CMd, control=mcprofileControl(), margin=NULL){
  if (is.null(rownames(CMn))) rownames(CMn) <- paste("C",1:nrow(CMn), sep="")
  if (is.null(colnames(CMn))) colnames(CMn) <- names(coefficients(object))
  if (ncol(CMn) != length(coefficients(object)) | ncol(CMd) != length(coefficients(object))) stop("Equal number of contrast and model coefficients needed!")
  if (control$fixed.range & length(margin) < 2) stop("2 margins have to be supplied to fix the step sizes between them (fixed.range = TRUE)!")
  kbbr <- function(x, Kn,Kd, b2, w){
    bw <- (b2*Kd[-w] %*% x - Kn[-w] %*% x) / (Kn[w] - b2*Kd[w])
    b <- numeric(length=length(Kn))
    b[(1:length(Kn))[-w]] <- x
    b[w] <- bw
    b
  }
  cutoff <- qnorm(1 - control$alphamax/nrow(CMn), 1)
  delta <- cutoff/control$steps
  prof <- stats:::profiler.nls(object)
  pars <- est <- prof$getFittedPars()
  estimate <- as.vector(CMn %*% est / CMd %*% est)
  SRDP <- list()
  for (i in 1:nrow(CMn)){
    Kn <- CMn[i,]
    Kd <- CMd[i,]
    neg <- Kn == 0 & Kd == 0
    w <- which(Kn != 0)[1]
    bstart <- est[-w]
    b <- b1 <- (Kn %*% est / Kd %*% est)[1,1]     
    bstart <- est[!neg][-w]
    par <- which(!neg)
    pars <- prof$getFittedPars()
    prof$setDefault(varying = par)
    mst <- 1
    while (z2 < cutoff & mst <= control$maxsteps){
      optfunc <- function(x){
        bw <- (b2*Kd[-w] %*% x - Kn[-w] %*% x) / (Kn[w] - b2*Kd[w])
        bopt <- numeric(length=length(Kn))
        bopt[(1:length(Kn))[-w]] <- x
        bopt[w] <- bw
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
      opar <- kbbr(conopt$par, Kn,Kd, b2, w)
      bstart <- opar[-w]
      z2 <- sqrt(conopt$value)
      z <- c(z, z2)
      b <- c(b, b2)
      stepsize <- sign(b2-b1)*(delta*((b2-b1)/(z2-z1)))
      b1 <- b2
      b2 <- b2 + stepsize
      z1 <- z2
      mst <- mst + 1
    }
    b1 <- (Kn %*% est / Kd %*% est)[1,1] 
    if (b1 != 0) b2 <- b1-abs(b1)*0.01 else b2 <- -0.01
    z1 <- z2 <- 0
    bstart <- est[!neg][-w]
    mst <- 1
    while (z2 < cutoff & mst <= control$maxsteps){
      optfunc <- function(x){
        bw <- (b2*Kd[-w] %*% x - Kn[-w] %*% x) / (Kn[w] - b2*Kd[w])
        bopt <- numeric(length=length(Kn))
        bopt[(1:length(Kn))[-w]] <- x
        bopt[w] <- bw
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
      opar <- kbbr(conopt$par, Kn,Kd, b2, w)
      bstart <- opar[-w]
      z2 <- sqrt(conopt$value)
      z <- c(z, -z2)
      b <- c(b, b2)
      stepsize <- sign(b2-b1)*(delta*((b2-b1)/(z2-z1)))
      b1 <- b2
      b2 <- b2 - stepsize
      z1 <- z2
      mst <- mst + 1
    }
    zd <- data.frame(z, b)
    mat <- zd[order(zd[,2]),]
    mat <- unique(na.omit(mat[!is.nan(mat[,1]),]))
    SRDP[[i]] <- mat
  }
  names(estimate) <-  names(SRDP) <- rownames(CMn)
  ########
  fsplist <- lapply(SRDP, function(z){
    try(interpSpline(z[,2], z[,1]), silent=TRUE)
  })
  bsplist <- lapply(SRDP, function(z){
    try(interpSpline(z[,1], z[,2]), silent=TRUE)
  })
  new(Class="mcprofileRatio", CMn=CMn, CMd=CMd, estimate=estimate, model=list(df=df.residual(object)), SRDP=SRDP, fsplines=fsplist, bsplines=bsplist, control=control, method="BFGS")
}
  


