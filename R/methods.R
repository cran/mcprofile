## Quadratic approximation
setGeneric("wald",function(object){standardGeneric("wald")})
setMethod(f="wald", signature="mcprofile", definition=function(object){
  SRDP <- object@SRDP
  est <- object@estimate
  cvu <- object@vest
  wsrdp <- lapply(1:length(SRDP), function(i){
    srdpi <- SRDP[[i]]
    b <- srdpi[,2]
    srdpi[,1] <- (b-est[i])/sqrt(cvu[i])
    srdpi 
  })
  names(wsrdp) <- names(SRDP)
  fsplist <- lapply(wsrdp, function(z){
    try(interpSpline(z[,2], z[,1]), silent=TRUE)
  })
  bsplist <- lapply(wsrdp, function(z){
    try(interpSpline(z[,1], z[,2]), silent=TRUE)
  })
  object@SRDP <- wsrdp
  object@fsplines <- fsplist
  object@bsplines <- bsplist
  return(object)
})

### higher order approximations
setGeneric("hoa",function(object){standardGeneric("hoa")})
setMethod(f="hoa", signature="mcprofile", definition=function(object){
  SRDP <- object@SRDP
  est <- object@estimate
  cvu <- object@vestunsc
  dv <- object@dvest
  disp <- object@model$dispersion
  hsrdp <- lapply(1:length(SRDP), function(i){
    srdpi <- SRDP[[i]]
    b <- srdpi[,2]
    r <- srdpi[,1]*sqrt(disp)
    q <- (b-est[i])/sqrt(cvu[i])
    j <- 1/dv[i]
    jp <- 1/cvu[i]
    j1 <- j/jp
    j0 <- 1/srdpi[,3]
    rho <- sqrt(j1/j0)
    srdpi[,1] <- (r + log(rho*q/r) / r)*1/sqrt(disp)
    srdpi <- srdpi[!is.nan(srdpi[,1]),]
    srdpi
  })
  monsrdp <- lapply(1:length(hsrdp), function(i){
    srdat <- hsrdp[[i]]
    ml <- srdat[srdat$b < est[i],]
    mll <- ml[1:(nrow(ml)-2),]
    mu <- srdat[srdat$b > est[i],]
    muu <- mu[3:nrow(mu),]
    rbind(mll, muu)
  })
  names(monsrdp) <- names(hsrdp) <- names(SRDP)  
  fsplist <- lapply(monsrdp, function(z){
    try(interpSpline(z[,2], z[,1]), silent=TRUE)
  })
  bsplist <- lapply(monsrdp, function(z){
    try(interpSpline(z[,1], z[,2]), silent=TRUE)
  })
  object@estimate <- unlist(lapply(bsplist, function(spl){
    ep <- try(predict(spl, 0)$y, silent=TRUE)
    if (class(ep) == "try-error") NA else ep
  }))
  object@SRDP <- hsrdp
  object@fsplines <- fsplist
  object@bsplines <- bsplist
  return(object)
})


#### confidence intervals
setMethod("confint", signature="mcprofile", definition=function(object, parm="missing", level=0.95, adjust="single-step", alternative="two.sided", quant=NULL){
  if (is.null(quant)){
    pam <- c("bonferroni", "none", "single-step")
    if (!(adjust %in% pam)) stop(paste("adjust has to be one of:", paste(pam, collapse=", ")))
    CM <- object@CM
    df <- object@model$df
    df <- if (!is.na(df)) df else NULL
    cr <- NULL
    if (adjust == "none" | nrow(CM) == 1){
      if (alternative == "two.sided"){
        alpha <- (1-level)/2
      } else {
        alpha <- 1-level
      }
      if (is.null(df)){
        quant <- qnorm(1-alpha)
      } else {
        quant <- qt(1-alpha, df=df)
      }
    }
    if (adjust == "bonferroni"){
      if (alternative == "two.sided"){
        alpha <- (1-level)/2
      } else {
        alpha <- 1-level
      }
      if (is.null(df)){
        quant <- qnorm(1-alpha/nrow(CM))
      } else {
        quant <- qt(1-alpha/nrow(CM), df=df)
      }
    }
    if (adjust == "single-step" & nrow(CM) > 1){
      require(mvtnorm)
      vc <- object@model$vcov
      VC <- CM %*% vc %*% t(CM)
      d <- 1/sqrt(diag(VC))
      dd <- diag(d)
      cr <- dd %*% VC %*% dd
      if (alternative == "two.sided"){
        if (is.null(df)){
          quant <- qmvnorm(level, corr=cr, tail="both.tails")$quantile
        } else {
          quant <- qmvt(level, df=df, corr=cr, tail="both.tails")$quantile
        }
      } else {
        if (is.null(df)){
          quant <- qmvnorm(level, corr=cr, tail="lower.tail")$quantile
        } else {
          quant <- qmvt(level, df=df, corr=cr, tail="lower.tail")$quantile
        }
      }
    }
  } else {
    adjust <- "user-defined"
  }
  if (alternative == "two.sided"){
    ci <- data.frame(t(sapply(object@bsplines, function(x, quant){
      pc <- try(c(predict(x, -quant)$y, predict(x, quant)$y), silent=TRUE)
      if (class(pc) == "try-error") c(NA, NA) else pc
    }, quant=quant)))
    names(ci) <- c("lower", "upper")
  }
  if (alternative == "less"){
    ci <- data.frame(sapply(object@bsplines, function(x, quant){
      pc <- try(rbind(predict(x, quant)$y), silent=TRUE)
      if (class(pc) == "try-error") rbind(c(NA)) else pc
    }, quant=quant))
    names(ci) <- "upper"
  }
  if (alternative == "greater"){
    ci <- data.frame(sapply(object@bsplines, function(x, quant){
      pc <- try(rbind(predict(x, -quant)$y), silent=TRUE)
      if (class(pc) == "try-error") rbind(c(NA)) else pc
    }, quant=quant))
    names(ci) <- "lower"
  }  
  new(Class="mcpconfint", confint=ci, quantile=quant, estimate=object@estimate, alternative=alternative, adjust=adjust)
})


setMethod("confint", signature="mcprofileRatio", definition=function(object, parm="missing", level=0.95, adjust="single-step", alternative="two.sided", quant=NULL){
  if (is.null(quant)){
    pam <- c("bonferroni", "none", "single-step")
    if (!(adjust %in% pam)) stop(paste("adjust has to be one of:", paste(pam, collapse=", ")))
    CMn <- object@CMn
    CMd <- object@CMd
    df <- object@model$df
    df <- if (!is.na(df)) df else NULL
    cr <- NULL
    if (adjust == "none" | nrow(CMn) == 1){
      if (alternative == "two.sided"){
        alpha <- (1-level)/2
      } else {
        alpha <- 1-level
      }
      if (is.null(df)){
        quant <- qnorm(1-alpha)
      } else {
        quant <- qt(1-alpha, df=df)
      }
    }
    if (adjust == "bonferroni"){
      if (alternative == "two.sided"){
        alpha <- (1-level)/2
      } else {
        alpha <- 1-level
      }
      if (is.null(df)){
        quant <- qnorm(1-alpha/nrow(CMn))
      } else {
        quant <- qt(1-alpha/nrow(CMn), df=df)
      }
    }
    if (adjust == "single-step" & nrow(CMn) > 1){
      require(mvtnorm)
      est <- object@estimate
      vc <- object@model$vcov
      cr <- matrix(rep(NA, nrow(CMn) * nrow(CMn)), nrow = nrow(CMn))
      for (i in 1:nrow(CMn)) {
        for (j in 1:nrow(CMn)) {
          cr[i,j] <- (est[i] * CMd[i,] - CMn[i,]) %*% vc %*% (est[j] * CMd[j,] - CMn[j,])/(sqrt((est[i] * CMd[i,] - CMn[i,]) %*% vc %*% (est[i] * CMd[i,] - CMn[i,])) * sqrt((est[j] * CMd[j,] - CMn[j,]) %*% vc %*% (est[j] * CMd[j,] - CMn[j,])))
        }
      }
      if (alternative == "two.sided"){
        if (is.null(df)){
          quant <- qmvnorm(level, corr=cr, tail="both.tails")$quantile
        } else {
          quant <- qmvt(level, df=df, corr=cr, tail="both.tails")$quantile
        }
      } else {
        if (is.null(df)){
          quant <- qmvnorm(level, corr=cr, tail="lower.tail")$quantile
        } else {
          quant <- qmvt(level, df=df, corr=cr, tail="lower.tail")$quantile
        }
      }
    }
  } else {
    adjust <- "user-defined"
  }
  if (alternative == "two.sided"){
    ci <- data.frame(t(sapply(object@bsplines, function(x, quant){
      pc <- try(c(predict(x, -quant)$y, predict(x, quant)$y), silent=TRUE)
      if (class(pc) == "try-error") c(NA, NA) else pc
    }, quant=quant)))
    names(ci) <- c("lower", "upper")
  }
  if (alternative == "less"){
    ci <- data.frame(sapply(object@bsplines, function(x, quant){
      pc <- try(rbind(predict(x, quant)$y), silent=TRUE)
      if (class(pc) == "try-error") rbind(c(NA)) else pc
    }, quant=quant))
    names(ci) <- "upper"
  }
  if (alternative == "greater"){
    ci <- data.frame(sapply(object@bsplines, function(x, quant){
      pc <- try(rbind(predict(x, -quant)$y), silent=TRUE)
      if (class(pc) == "try-error") rbind(c(NA)) else pc
    }, quant=quant))
    names(ci) <- "lower"
  }  
  new(Class="mcpconfint", confint=ci, quantile=quant, estimate=object@estimate, alternative=alternative, adjust=adjust)
})






# exp ci
setMethod(f="exp", signature="mcpconfint", definition=function(x){
  x@estimate <- exp(x@estimate)
  names(x@estimate) <- paste("exp( ", names(x@estimate), " )", sep="")
  x@confint <- exp(x@confint)
  return(x)
})

# expit ci
setMethod(f="expit", signature="mcpconfint", definition=function(x){
  x@estimate <- exp(x@estimate)/(1 + exp(x@estimate))
  names(x@estimate) <- paste("expit( ", names(x@estimate), " )", sep="")
  x@confint <- exp(x@confint)/(1 + exp(x@confint))
  return(x)
})



## Tests
setGeneric("test",function(object, adjust="bonferroni", alternative="two.sided", margin=0){standardGeneric("test")})
setMethod("test", signature="mcprofile", definition=function(object, adjust="single-step", alternative="two.sided", margin=0){
  CM <- object@CM
  est <- object@estimate
  df <- object@model$df
  df <- if (!is.na(df)) df else NULL
  if (!(alternative %in% c("two.sided", "less", "greater"))) stop("alternative has to be one of 'two.sided', 'less', or 'greater'")
  pam <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")  
  if (!(adjust %in% c(pam, "single-step"))) stop(paste("adjust has to be one of:", paste(c(pam, "single-step"), collapse=", ")))
  ##### Find stats
  ptest <- function(ispl, delta){
    pst <- try(predict(ispl, delta)$y, silent=TRUE)
    if (class(pst) == "try-error") NA else pst   
  }
  stat <- sapply(object@fsplines, function(x) ptest(x, delta=margin))
  naid <- is.na(stat)
  stat[is.na(stat)] <- 0
  if (is.null(df)){
    switch(alternative, less = {
      praw <- pnorm(stat, lower.tail=TRUE)
    }, greater = {
      praw <- pnorm(stat, lower.tail=FALSE)
    }, two.sided = {
      praw <- pmin(1, pnorm(abs(stat), lower.tail=FALSE)*2)
    })
  } else {
    switch(alternative, less = {
      praw <- pt(stat, df=df, lower.tail=TRUE)
    }, greater = {
      praw <- pt(stat, df=df, lower.tail=FALSE)
    }, two.sided = {
      praw <- pmin(1, pt(abs(stat), df=df, lower.tail=FALSE)*2)
    })
  }
  if (length(praw) > 1){
    if (adjust %in% pam) padj <- p.adjust(praw, method=adjust)
    pfct <- function(q) {
      switch(alternative, two.sided = {
      	low <- rep(-abs(q), dim)
      	upp <- rep(abs(q), dim)
      }, less = {
      	low <- rep(q, dim)
      	upp <- rep(Inf, dim)
      }, greater = {
      	low <- rep(-Inf, dim)
      	upp <- rep(q, dim)
      })
      if (is.null(df)){
      	pmvnorm(lower = low, upper = upp, corr = cr)
      } else {
      	pmvt(lower = low, upper = upp, df=df, corr = cr)
      }
    }
    if (adjust == "single-step"){
      require(mvtnorm)
      vc <- object@model$vcov
      VC <- CM %*% vc %*% t(CM)
      d <- 1/sqrt(diag(VC))
      dd <- diag(d)
      cr <- dd %*% VC %*% dd
      dim <- ncol(cr)
      padj <- numeric(length(stat))
      error <- 0
      for (i in 1:length(stat)) {
      	tmp <- pfct(stat[i])
      	if (error < attr(tmp, "error")) error <- attr(tmp, "error")
      	padj[i] <- tmp
      }
      padj <- 1 - padj
      attr(padj, "error") <- error
    }
  } else { padj <- praw }
  stat[naid] <- NA
  padj[naid] <- NA
  new(Class="mcptest", stat = stat, pvalues = padj, margin = margin, estimate = est, alternative = alternative, adjust=adjust)  
})


setMethod("test", signature="mcprofileRatio", definition=function(object, adjust="single-step", alternative="two.sided", margin=1){
  CMn <- object@CMn
  CMd <- object@CMd
  est <- object@estimate
  df <- object@model$df
  df <- if (!is.na(df)) df else NULL
  if (!(alternative %in% c("two.sided", "less", "greater"))) stop("alternative has to be one of 'two.sided', 'less', or 'greater'")
  pam <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")  
  if (!(adjust %in% c(pam, "single-step"))) stop(paste("adjust has to be one of:", paste(c(pam, "single-step"), collapse=", ")))
  ##### Find stats
  ptest <- function(ispl, delta){
    pst <- try(predict(ispl, delta)$y, silent=TRUE)
    if (class(pst) == "try-error") NA else pst   
  }
  stat <- sapply(object@fsplines, function(x) ptest(x, delta=margin))
  naid <- is.na(stat)
  stat[is.na(stat)] <- 0
  if (is.null(df)){
    switch(alternative, less = {
      praw <- pnorm(stat, lower.tail=TRUE)
    }, greater = {
      praw <- pnorm(stat, lower.tail=FALSE)
    }, two.sided = {
      praw <- pmin(1, pnorm(abs(stat), lower.tail=FALSE)*2)
    })
  } else {
    switch(alternative, less = {
      praw <- pt(stat, df=df, lower.tail=TRUE)
    }, greater = {
      praw <- pt(stat, df=df, lower.tail=FALSE)
    }, two.sided = {
      praw <- pmin(1, pt(abs(stat), df=df, lower.tail=FALSE)*2)
    })
  }
  if (length(praw) > 1){
    if (adjust %in% pam) padj <- p.adjust(praw, method=adjust)
    pfct <- function(q) {
      switch(alternative, two.sided = {
      	low <- rep(-abs(q), dim)
      	upp <- rep(abs(q), dim)
      }, less = {
      	low <- rep(q, dim)
      	upp <- rep(Inf, dim)
      }, greater = {
      	low <- rep(-Inf, dim)
      	upp <- rep(q, dim)
      })
      if (is.null(df)){
      	pmvnorm(lower = low, upper = upp, corr = cr)
      } else {
      	pmvt(lower = low, upper = upp, df=df, corr = cr)
      }
    }
    if (adjust == "single-step"){
      require(mvtnorm)
      est <- object@estimate
      vc <- object@model$vcov
      cr <- matrix(rep(NA, nrow(CMn) * nrow(CMn)), nrow = nrow(CMn))
      for (i in 1:nrow(CMn)) {
        for (j in 1:nrow(CMn)) {
          cr[i,j] <- (est[i] * CMd[i,] - CMn[i,]) %*% vc %*% (est[j] * CMd[j,] - CMn[j,])/(sqrt((est[i] * CMd[i,] - CMn[i,]) %*% vc %*% (est[i] * CMd[i,] - CMn[i,])) * sqrt((est[j] * CMd[j,] - CMn[j,]) %*% vc %*% (est[j] * CMd[j,] - CMn[j,])))
        }
      }
      dim <- ncol(cr)
      padj <- numeric(length(stat))
      error <- 0
      for (i in 1:length(stat)) {
      	tmp <- pfct(stat[i])
      	if (error < attr(tmp, "error")) error <- attr(tmp, "error")
      	padj[i] <- tmp
      }
      padj <- 1 - padj
      attr(padj, "error") <- error
    }
  } else { padj <- praw }
  stat[naid] <- NA
  padj[naid] <- NA
  new(Class="mcptest", stat = stat, pvalues = padj, margin = margin, estimate = est, alternative = alternative, adjust=adjust)  
})

