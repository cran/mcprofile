##################################
mcpcalcmer <- function(object, CM, control=mcprofileControl(), margin=NULL){
  #######
  glmer_finalizeDG <- function(fr, FL, glmFit, start, nAGQ, verbose){
    if (is.list(start) && all(sort(names(start)) == sort(names(FL)))) start <- list(ST = start)
    if (is.numeric(start)) start <- list(STpars = start)
    dm <- lme4:::mkZt(FL, start[["ST"]])
    ft <- lme4:::famType(glmFit$family)
    dm$dd[names(ft)] <- ft
    useSc <- as.integer(!(lme4:::famNms[dm$dd[["fTyp"]]] %in% c("binomial", "poisson")))
    dm$dd[["useSc"]] <- useSc
    M1 <- length(levels(dm$flist[[1]]))
    n <- ncol(dm$Zt)
    if ((nAGQ <- as.integer(nAGQ)) < 1) nAGQ <- 1L
    if (nAGQ%%2 == 0) nAGQ <- nAGQ + 1L
    dm$dd["nAGQ"] <- as.integer(nAGQ)
    AGQlist <- .Call(lme4:::lme4_ghq, nAGQ)
    y <- unname(as.double(glmFit$y))
    p <- dm$dd[["p"]]
    dm$dd["verb"] <- as.integer(verbose)
    fixef <- fr$fixef
    fixef[] <- coef(glmFit)
    if (!is.null(ff <- start$fixef) && is.numeric(ff) && length(ff) == length(fixef)) fixef <- ff
    ans <- new(Class = "mer", env = new.env(), nlmodel = (~I(x))[[2]], frame = fr$mf, call = call("foo"), flist = dm$flist, Zt = dm$Zt, X = fr$X, y = y, pWt = unname(glmFit$prior.weights), offset = unname(fr$off), Gp = unname(dm$Gp), dims = dm$dd, ST = dm$ST, A = dm$A, Cm = dm$Cm, Cx = (dm$A)@x, L = dm$L, deviance = dm$dev, fixef = fixef, ranef = numeric(dm$dd[["q"]]), u = numeric(dm$dd[["q"]]), eta = unname(glmFit$linear.predictors), mu = unname(glmFit$fitted.values), muEta = numeric(dm$dd[["n"]]), var = numeric(dm$dd[["n"]]), resid = unname(glmFit$residuals), sqrtXWt = as.matrix(numeric(dm$dd[["n"]])), sqrtrWt = numeric(dm$dd[["n"]]), RZX = matrix(0, dm$dd[["q"]], p), RX = matrix(0, p, p), ghx = AGQlist[[1]], ghw = AGQlist[[2]])
    if (!is.null(stp <- start$STpars) && is.numeric(stp)) {
      STp <- .Call(lme4:::mer_ST_getPars, ans)
      if (length(STp) == length(stp))
        .Call(lme4:::mer_ST_setPars, ans, stp)
    }
    lme4:::mer_finalize(ans)
    ans
  }
  ########
  if (is.null(rownames(CM))) rownames(CM) <- paste("C",1:nrow(CM), sep="")
  if (is.null(colnames(CM))) colnames(CM) <- names(fixef(object))
  if (ncol(CM) != length(fixef(object))) stop("Equal number of contrast and model coefficients needed!")
  if (control$fixed.range & length(margin) < 2) stop("2 margins have to be supplied to fix the step sizes between them (fixed.range = TRUE)!")
  mc <- object@call
  family <- eval(mc$family)
  if (is.null(family)) family <- gaussian()
  form <- mc$formula
  gauss <- FALSE
  if (family$family == "gaussian" && family$link == "identity") {
    mc[[1]] <- as.name("lmer")
    mc$family <- NULL
    gauss <- TRUE
    #return(eval.parent(mc))
  }

  fr <- lme4:::lmerFrames(mc, form, NULL)
  offset <- wts <- NULL
  if (length(fr$wts)) wts <- fr$wts else wts <- NULL
  if (length(fr$off)) offset <- fr$off else offset <- NULL
  cv <- do.call(lme4:::lmerControl, list())
  verbose <- cv$msVerbose

  dev0 <- object@deviance[1]
  mod <- list(coefficients = fixef(object), vcov = as.matrix(vcov(object)), X = object@X, deviance = dev0, response = fr$Y, weights = wts, family = family, fitted = fitted(object), offset = offset, df = NA, converged = TRUE)
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
  off <- fr$off
  X <- fr$X
  y <- fr$Y
  est <- fixef(object)
  estimate <- as.vector(CM %*% est)
  vest <- diag(CM %*% as.matrix(vcov(object)) %*% t(CM))
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
    fr$X <- X[,neg, drop=FALSE]
    fr$fixef <- est[neg]
    bp <- est[pos]
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
        fr$off <- as.numeric(X[,!neg,drop=FALSE] %*% bopt[!neg])
        startc <- get("startc", envir=.start)
        FL <- lme4:::lmerFactorList(form, fr, 0L, 0L)
        FL$dims["mxit"] <- cv$maxIter
        FL$dims["mxfn"] <- cv$maxFN
        if (gauss == TRUE){
          ans <- list(fr = fr, FL = FL, start = NULL, REML = FALSE, verbose = FALSE)
          ans <- do.call(lme4:::lmer_finalize, ans)
        } else {
          glmFit <- glm.fit(fr$X, fr$Y, weights = wts, offset = fr$off, family = family, intercept = attr(attr(fr$mf, "terms"), "intercept") > 0)
          nAGQ <- 1
          ans <- list(fr = fr, FL = FL, glmFit = glmFit, start = NULL, nAGQ = nAGQ, verbose = verbose)
          ans <- do.call(glmer_finalizeDG, ans)
        }
        ans@call <- mc
        #assign("startc", cfm$coefficients, envir = .start)
        ans@deviance[1]
      }
      conopt <- try(optim(bstart, optd, method="BFGS"), silent=TRUE)
      if (class(conopt) == "try-error") break
      opar <- kbb(conopt$par, K, b2, w)
      bstart <- opar[-w]
      dev1 <- conopt$value
      z2 <- sqrt(abs(dev1 - dev0))
      z <- c(z, z2)
      b <- c(b, b2)
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
        fr$off <- as.numeric(X[,!neg,drop=FALSE] %*% bopt[!neg])
        FL <- lme4:::lmerFactorList(form, fr, 0L, 0L)
        FL$dims["mxit"] <- cv$maxIter
        FL$dims["mxfn"] <- cv$maxFN
        if (gauss == TRUE){
          ans <- list(fr = fr, FL = FL, start = NULL, REML = FALSE, verbose = FALSE)
          ans <- do.call(lme4:::lmer_finalize, ans)
        } else {
          glmFit <- glm.fit(fr$X, fr$Y, weights = wts, offset = fr$off, family = family, intercept = attr(attr(fr$mf, "terms"), "intercept") > 0)
          nAGQ <- 1
          ans <- list(fr = fr, FL = FL, glmFit = glmFit, start = NULL, nAGQ = nAGQ, verbose = verbose)
          ans <- do.call(glmer_finalizeDG, ans)
        }
        ans@call <- mc
        #assign("startc", cfm$coefficients, envir = .start)
        ans@deviance[1]
      }
      conopt <- try(optim(bstart, optd, method="BFGS"), silent=TRUE)
      if (class(conopt) == "try-error") break
      opar <- kbb(conopt$par, K, b2, w)
      bstart <- opar[-w]      
      dev1 <- conopt$value
      z2 <- sqrt(abs(dev1 - dev0))
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
    Kv <- rbind(K, nuiKv[apply(nuiKv,1,sum) != 0,])
  }
  names(SRDP) <- rownames(CM)
  ########
  fsplist <- lapply(SRDP, function(z){
    try(interpSpline(z[,2], z[,1]), silent=TRUE)
  })
  new(Class="mcprofile", CM=CM, estimate=estimate, vest=vest, vestunsc=numeric(), dvest=numeric(), model=mod, SRDP=SRDP, fsplines=fsplist, control=control, method="BFGS")
}


######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################




mcpcalcmerlmer <- function(object, CM, control=mcprofileControl(), margin=NULL){
  glmer_finalizeDG <- function(fr, FL, glmFit, start, nAGQ, verbose){
    if (is.list(start) && all(sort(names(start)) == sort(names(FL)))) start <- list(ST = start)
    if (is.numeric(start)) start <- list(STpars = start)
    dm <- lme4:::mkZt(FL, start[["ST"]])
    ft <- lme4:::famType(glmFit$family)
    dm$dd[names(ft)] <- ft
    useSc <- as.integer(!(lme4:::famNms[dm$dd[["fTyp"]]] %in% c("binomial", "poisson")))
    dm$dd[["useSc"]] <- useSc
    M1 <- length(levels(dm$flist[[1]]))
    n <- ncol(dm$Zt)
    if ((nAGQ <- as.integer(nAGQ)) < 1) nAGQ <- 1L
    if (nAGQ%%2 == 0) nAGQ <- nAGQ + 1L
    dm$dd["nAGQ"] <- as.integer(nAGQ)
    AGQlist <- .Call(lme4:::lme4_ghq, nAGQ)
    y <- unname(as.double(glmFit$y))
    p <- dm$dd[["p"]]
    dm$dd["verb"] <- as.integer(verbose)
    fixef <- fr$fixef
    fixef[] <- coef(glmFit)
    if (!is.null(ff <- start$fixef) && is.numeric(ff) && length(ff) == length(fixef)) fixef <- ff
    ans <- new(Class = "mer", env = new.env(), nlmodel = (~I(x))[[2]], frame = fr$mf, call = call("foo"), flist = dm$flist, Zt = dm$Zt, X = fr$X, y = y, pWt = unname(glmFit$prior.weights), offset = unname(fr$off), Gp = unname(dm$Gp), dims = dm$dd, ST = dm$ST, A = dm$A, Cm = dm$Cm, Cx = (dm$A)@x, L = dm$L, deviance = dm$dev, fixef = fixef, ranef = numeric(dm$dd[["q"]]), u = numeric(dm$dd[["q"]]), eta = unname(glmFit$linear.predictors), mu = unname(glmFit$fitted.values), muEta = numeric(dm$dd[["n"]]), var = numeric(dm$dd[["n"]]), resid = unname(glmFit$residuals), sqrtXWt = as.matrix(numeric(dm$dd[["n"]])), sqrtrWt = numeric(dm$dd[["n"]]), RZX = matrix(0, dm$dd[["q"]], p), RX = matrix(0, p, p), ghx = AGQlist[[1]], ghw = AGQlist[[2]])
    if (!is.null(stp <- start$STpars) && is.numeric(stp)) {
      STp <- .Call(lme4:::mer_ST_getPars, ans)
      if (length(STp) == length(stp))
        .Call(lme4:::mer_ST_setPars, ans, stp)
    }
    lme4:::mer_finalize(ans)
    ans
  }
  if (is.null(rownames(CM))) rownames(CM) <- paste("C",1:nrow(CM), sep="")
  if (is.null(colnames(CM))) colnames(CM) <- names(fixef(object))
  if (ncol(CM) != length(fixef(object))) stop("Equal number of contrast and model coefficients needed!")
  if (control$fixed.range & length(margin) < 2) stop("2 margins have to be supplied to fix the step sizes between them (fixed.range = TRUE)!")
  mc <- object@call
  nAGQ <- mc$nAGQ
  if (is.null(nAGQ)) nAGQ <- 1
  family <- eval(mc$family)
  if (is.null(family)) family <- gaussian()
  form <- mc$formula
  gauss <- FALSE
  if (family$family == "gaussian" && family$link == "identity") {
    mc[[1]] <- as.name("lmer")
    mc$family <- NULL
    gauss <- TRUE
    #return(eval.parent(mc))
  }

  fr <- lme4:::lmerFrames(mc, form, NULL)
  offset <- wts <- NULL
  if (length(fr$wts)) wts <- fr$wts else wts <- NULL
  if (length(fr$off)) offset <- fr$off else offset <- NULL
  cv <- do.call(lme4:::lmerControl, list())
  verbose <- cv$msVerbose

  mod <- list(coefficients = fixef(object), vcov = as.matrix(vcov(object)), X = object@X, deviance = object@deviance[1], response = fr$Y, weights = wts, family = family, fitted = fitted(object), offset = offset,   df = NA, converged = TRUE)
  kbb <- function(x, k, p, w){
    K <- -1*c(-p, k[-w])/k[w]
    hv <- c(1,x)
    b1 <- (K %*% hv)[1,1]
    b <- numeric(length=length(k))
    b[(1:length(k))[-w]] <- x
    b[w] <- b1
    b
  }
  lmeroff <- function(fr, off, gauss, start, family, nAGQ, startf){
    fr$fixef <- startf
    fr$off <- off
    cv <- do.call(lme4:::lmerControl, list())
    verbose <- cv$msVerbose
    FL <- lme4:::lmerFactorList(form, fr, 0L, 0L)
    FL$dims["mxit"] <- cv$maxIter
    FL$dims["mxfn"] <- cv$maxFN
    if (gauss == TRUE){
      ans <- list(fr = fr, FL = FL, start = start, REML = FALSE, verbose = FALSE)
      ans <- do.call(lme4:::lmer_finalize, ans)
    } else {
      glmFit <- glm.fit(fr$X, fr$Y, weights = NULL, offset = off, family = family, intercept = attr(attr(fr$mf, "terms"), "intercept") > 0)
      ans <- list(fr = fr, FL = FL, glmFit = glmFit, start = start, nAGQ = nAGQ, verbose = verbose)
      ans <- do.call(glmer_finalizeDG, ans)
    }
    ans@call <- mc
    ans
  }
  cutoff <- qnorm(1 - control$alphamax/nrow(CM), 1)
  delta <- cutoff/control$steps
  off <- fr$off
  X <- fr$X
  y <- fr$Y
  dev0 <- object@deviance[1]
  est <- fixef(object)
  .start <- new.env()
  SRDP <- list()
  estimate <- numeric(length=nrow(CM))
  vest <- numeric(length=nrow(CM))
  for (i in 1:nrow(CM)){
    K <- CM[i,]
    if (any(K < 0) & any(K > 0)){
      Xi <- makeDesign(K, mod$X)
      startf <- c(mod$coefficients[1], K %*% mod$coefficients, mod$coefficients[K == 0])
    }
    if (all(K >= 0)){
      Xi <- pmakeDesign(K, mod$X)
      startf <- c(mod$coefficients[K == 0][1], K %*% mod$coefficients, mod$coefficients[K == 0][-1])
    }
    Xio <- Xi[,-2,drop=FALSE]
    Xib <- Xi[,2,drop=FALSE]
    fr$X <- Xi
    fr$fixef <- startf
    FL <- lme4:::lmerFactorList(form, fr, 0L, 0L)
    FL$dims["mxit"] <- cv$maxIter
    FL$dims["mxfn"] <- cv$maxFN
    if (gauss == TRUE){
      ans <- list(fr = fr, FL = FL, start = NULL, REML = FALSE, verbose = FALSE)
      ans <- do.call(lme4:::lmer_finalize, ans)
    } else {
      glmFit <- glm.fit(fr$X, fr$Y, weights = wts, offset = NULL, family = family, intercept = attr(attr(fr$mf, "terms"), "intercept") > 0)
      nAGQ <- 1
      ans <- list(fr = fr, FL = FL, glmFit = glmFit, start = unlist(object@ST), nAGQ = nAGQ, verbose = verbose)
      ans <- do.call(glmer_finalizeDG, ans)
    }
    ans@call <- mc
    pans <- oans <- ans
    ddev <- oans@deviance[1]
    estimate[i] <- pp1 <- oans@fixef[2]
    vest[i] <- as.matrix(vcov(oans))[2,2]
    fr$X <- Xio
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
    startc <- unlist(ans@ST)
    startf <- fixef(ans)[-2]
    mst <- 1
    while((zp1 < cutoff | pp1 < maxmargin) & mst <= control$maxsteps){
      pans <- lmeroff(fr, (Xib %*% pp2)[,1], start=startc, gauss=gauss, family, nAGQ, startf)
      zp <- sqrt(abs(pans@deviance[1] - ddev))
      if (!control$fixed.range) stepsize <- sign(pp2-pp1)*(delta*((pp2-pp1)/(zp-zp1)))
      pp <- pp2 + stepsize
      pp1 <- pp2
      pp2 <- pp
      zp1 <- zp
      z <- c(z, zp1)
      b <- c(b, pp1)
      startc <- unlist(pans@ST)
      startf <- fixef(pans)
      if (control$fixed.range & length(b) == control$steps+1) break
      mst <- mst + 1
    }
    pans <- oans
    pp1 <- oans@fixef[2]
    if (pp1 != 0) pp2 <- pp1-abs(pp1*0.01) else pp2 <- -0.01
    if (control$fixed.range){
      stepsize <- diff(c(min(margin), pp1)) / control$steps
      pp2 <- pp1 - stepsize
    }    
    zp1 <- 0
    startc <- unlist(ans@ST)
    startf <- fixef(ans)[-2]
    mst <- 1
    while((zp1 < cutoff | pp1 > maxmargin) & mst <= control$maxsteps){
      pans <- lmeroff(fr, (Xib %*% pp2)[,1], start=startc, gauss=gauss, family, nAGQ, startf)
      zp <- sqrt(abs(pans@deviance[1] - ddev))
      if (!control$fixed.range) stepsize <- sign(pp2-pp1)*(delta*((pp2-pp1)/(zp-zp1)))
      pp <- pp2 - stepsize
      pp1 <- pp2
      pp2 <- pp
      zp1 <- zp
      z <- c(z, -zp1)
      b <- c(b, pp1)
      startc <- unlist(pans@ST)
      startf <- fixef(pans)
      if (control$fixed.range & length(b) == 2*control$steps+1) break
      mst <- mst + 1
    }
 
    zd <- data.frame(z, b)
    mat <- zd[order(zd[,2]),]
    mat <- unique(na.omit(mat[!is.nan(mat[,1]),]))
    SRDP[[i]] <- mat
  }
  names(SRDP) <- rownames(CM)
  ########
  fsplist <- lapply(SRDP, function(z){
    try(interpSpline(z[,2], z[,1]), silent=TRUE)
  })
  new(Class="mcprofile", CM=CM, estimate=estimate, vest=vest, vestunsc=numeric(), dvest=numeric(), model=mod, SRDP=SRDP, fsplines=fsplist, control=control, method="lmer")
}

