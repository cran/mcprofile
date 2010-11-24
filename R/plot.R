setMethod("plot", signature=c("mcprofile"), definition=function(x, y, squared=FALSE, trans=identity, ...){
  dp <- x@SRDP
  fs <- x@fsplines
  z <- sapply(dp, function(x) x[,1])
  b <- sapply(dp, function(x) x[,2])
  d <- data.frame(b=as.vector(unlist(b)), z=as.vector(unlist(z)))
  ff <- as.factor(rep(names(dp), unlist(lapply(dp, nrow))))
  pspl <- lapply(1:length(fs), function(i){
    pseq <- seq(min(dp[[i]][,2]), max(dp[[i]][,2]), length=100)
    psp <- try(predict(fs[[i]], pseq), silent=TRUE)
    if (class(psp) == "try-error")  cbind(rep(NA, length(pseq)), rep(NA, length(pseq))) else cbind(psp$x, psp$y)
  })
  sx <- sapply(pspl, function(x) x[,1])
  sy <- sapply(pspl, function(x) x[,2])
  sdat <- data.frame(x=as.vector(sx), y=as.vector(sy))
  ffs <- as.factor(rep(names(dp), unlist(lapply(pspl, nrow))))
  ddd <- data.frame(x=c(d[,1],sdat[,1]), y=c(d[,2], sdat[,2]), grp=c(as.character(ff), as.character(ffs)), typ=c(rep("obs", nrow(d)), rep("preds", nrow(sdat))))
  if (squared == TRUE) ddd$y <- ddd$y^2
  ddd$x <- trans(ddd$x)
  pg <- function(x, y, group.number, ...){
    panel.abline(h=0, col="lightgrey")
    if (group.number == 1){
      panel.segments(x, 0, x, y, col="black")
    }
    if (group.number == 2){
      panel.lines(x, y, col="blue2", lwd=2)
    }    
  }
  xyplot(y ~ x | grp, group=typ, data=ddd, panel=panel.superpose, panel.groups=pg, scale=list(x=list(relation="free")), xlab=expression(beta), ylab="z", ...)
})

setMethod("plot", signature=c("mcprofileRatio"), definition=function(x, y, squared=FALSE, trans=identity, ...){
  dp <- x@SRDP
  fs <- x@fsplines
  z <- sapply(dp, function(x) x[,1])
  b <- sapply(dp, function(x) x[,2])
  d <- data.frame(b=as.vector(unlist(b)), z=as.vector(unlist(z)))
  ff <- as.factor(rep(names(dp), unlist(lapply(dp, nrow))))
  pspl <- lapply(1:length(fs), function(i){
    pseq <- seq(min(dp[[i]][,2]), max(dp[[i]][,2]), length=100)
    psp <- try(predict(fs[[i]], pseq), silent=TRUE)
    if (class(psp) == "try-error")  cbind(rep(NA, length(pseq)), rep(NA, length(pseq))) else cbind(psp$x, psp$y)
  })
  sx <- sapply(pspl, function(x) x[,1])
  sy <- sapply(pspl, function(x) x[,2])
  sdat <- data.frame(x=as.vector(sx), y=as.vector(sy))
  ffs <- as.factor(rep(names(dp), unlist(lapply(pspl, nrow))))
  ddd <- data.frame(x=c(d[,1],sdat[,1]), y=c(d[,2], sdat[,2]), grp=c(as.character(ff), as.character(ffs)), typ=c(rep("obs", nrow(d)), rep("preds", nrow(sdat))))
  if (squared == TRUE) ddd$y <- ddd$y^2
  ddd$x <- trans(ddd$x)
  pg <- function(x, y, group.number, ...){
    panel.abline(h=0, col="lightgrey")
    if (group.number == 1){
      panel.segments(x, 0, x, y, col="black")
    }
    if (group.number == 2){
      panel.lines(x, y, col="blue2", lwd=2)
    }    
  }
  xyplot(y ~ x | grp, group=typ, data=ddd, panel=panel.superpose, panel.groups=pg, scale=list(x=list(relation="free")), xlab=expression(beta), ylab="z", ...)
})




setMethod("plot", signature="mcpconfint", definition=function(x, y, margin=NULL, pch=19, trans=identity, ...){
  estimate <- trans(x@estimate)
  if (x@alternative == "two.sided"){
    lower <- trans(x@confint[,1])
    upper <- trans(x@confint[,2])
  }
  if (x@alternative == "greater"){
    lower <- trans(x@confint[,1])
    upper <- NULL
  }
  if (x@alternative == "less"){
    upper <- trans(x@confint[,1])
    lower <- NULL
  }  
  n <- length(estimate)
  enames <- names(estimate)
  alle <- na.omit(c(estimate, lower, upper, margin))
  alle <- alle[alle > -Inf & alle < Inf]
  if (length(alle) == 0) alle <- c(-1,1)
  mini <- min(alle)
  maxi <- max(alle)
  dis <- abs(maxi-mini)
  if (dis == 0) dis <- 1
  llim <- mini-0.05*dis
  ulim <- maxi+0.05*dis
  plot.new()
  par(usr=c(llim, ulim, 0.5, n+0.5)) 
  axis(2, at=n:1, labels=enames, las=2, pos=llim-0.01*dis, lwd=0)
  axis(1)
  for (i in n:1){
    if (is.null(lower)){
      ll <- llim
    } else {
      if (lower[i] == -Inf | is.na(lower[i])){
        ll <- llim
      } else {
        ll <- lower[i]
      }
    }
    if (is.null(upper)){
      ul <- ulim
    } else {
      if (upper[i] == Inf | is.na(upper[i])){
        ul <- ulim
      } else {
        ul <- upper[i]
      }
    }
    acode <- 3
    if (ll == llim) acode <- 2
    if (ul == ulim) acode <- 1
    if (ul == ulim & ll == llim) acode <- 0
    abline(h=n-i+1, col="grey", lty=3)
    arrows(ll,n-i+1,ul,n-i+1, code=acode, angle=90, length=0.1, ...)
    if (!is.na(estimate[i]) & estimate[i] > -Inf & estimate[i] < Inf) points(estimate[i],n-i+1, pch=pch, ...)
  }
})



setMethod("plot", signature=c(x="mcprofile", y="mcpconfint"), definition=function(x, y, squared=FALSE, trans=identity, ...){
  dp <- x@SRDP
  fs <- x@fsplines
  pspl <- lapply(1:length(fs), function(i){
    pseq <- seq(min(dp[[i]][,2]), max(dp[[i]][,2]), length=100)
    psp <- try(predict(fs[[i]], pseq), silent=TRUE)
    if (class(psp) == "try-error")  cbind(rep(NA, length(pseq)), rep(NA, length(pseq))) else cbind(psp$x, psp$y)
  })
  sx <- sapply(pspl, function(x) x[,1])
  sy <- sapply(pspl, function(x) x[,2])
  sdat <- data.frame(x=as.vector(sx), y=as.vector(sy))
  ffs <- as.factor(rep(names(dp), unlist(lapply(pspl, nrow))))

  cc <- unlist(y@confint)
  ffci <- as.factor(rep(names(dp), ncol(y@confint)))
  qul <- -y@quantile
  quu <- y@quantile  
  if (y@alternative == "two.sided") quant <- c(rep(-y@quantile, length(cc)/2), rep(y@quantile, length(cc)/2))
  if (y@alternative == "less"){
    quant <- rep(y@quantile, length(cc))
    qul <- NULL
  }
  if (y@alternative == "greater"){
    quant <- rep(-y@quantile, length(cc))
    quu <- NULL
  }
  est <- y@estimate
  zer <- rep(0, length(est))
  ffest <- as.factor(names(dp))
  
  ddd <- data.frame(x=c(cc, est ,sdat[,1]), y=c(quant, zer, sdat[,2]), grp=c(as.character(ffci), as.character(ffest), as.character(ffs)), typ=c(rep("ci", length(cc)), rep("est", length(est)), rep("preds", nrow(sdat))))
  if (squared == TRUE){
    ddd$y <- ddd$y^2
    qul <- quu <- y@quantile^2
  }
  ymin <- min(ddd[,2][!is.na(ddd[,2])])
  ddd$x <- trans(ddd$x)
  pg <- function(x, y, group.number, ...){
    panel.abline(h=0, col="lightgrey")
    if (!is.null(qul)) panel.abline(h=qul, col="grey", lty=2)
    if (!is.null(quu)) panel.abline(h=quu, col="grey", lty=2)     
    if (group.number == 1){
      panel.segments(x, ymin, x, y, col="black", lwd=2)
    }
    if (group.number == 2){
      panel.segments(x, ymin, x, y, col="black", lwd=2, lty=3)
    }    
    if (group.number == 3){
      panel.lines(x, y, col="blue2", lwd=2)
    }    
  }
  xyplot(y ~ x | grp, group=typ, data=ddd, panel=panel.superpose, panel.groups=pg, scale=list(x=list(relation="free")), xlab=expression(beta), ylab="z", ...)
})

setMethod("plot", signature=c(x="mcprofileRatio", y="mcpconfint"), definition=function(x, y, squared=FALSE, trans=identity, ...){
  dp <- x@SRDP
  fs <- x@fsplines
  pspl <- lapply(1:length(fs), function(i){
    pseq <- seq(min(dp[[i]][,2]), max(dp[[i]][,2]), length=100)
    psp <- try(predict(fs[[i]], pseq), silent=TRUE)
    if (class(psp) == "try-error")  cbind(rep(NA, length(pseq)), rep(NA, length(pseq))) else cbind(psp$x, psp$y)
  })
  sx <- sapply(pspl, function(x) x[,1])
  sy <- sapply(pspl, function(x) x[,2])
  sdat <- data.frame(x=as.vector(sx), y=as.vector(sy))
  ffs <- as.factor(rep(names(dp), unlist(lapply(pspl, nrow))))

  cc <- unlist(y@confint)
  ffci <- as.factor(rep(names(dp), ncol(y@confint)))
  qul <- -y@quantile
  quu <- y@quantile  
  if (y@alternative == "two.sided") quant <- c(rep(-y@quantile, length(cc)/2), rep(y@quantile, length(cc)/2))
  if (y@alternative == "less"){
    quant <- rep(y@quantile, length(cc))
    qul <- NULL
  }
  if (y@alternative == "greater"){
    quant <- rep(-y@quantile, length(cc))
    quu <- NULL
  }
  est <- y@estimate
  zer <- rep(0, length(est))
  ffest <- as.factor(names(dp))
  
  ddd <- data.frame(x=c(cc, est ,sdat[,1]), y=c(quant, zer, sdat[,2]), grp=c(as.character(ffci), as.character(ffest), as.character(ffs)), typ=c(rep("ci", length(cc)), rep("est", length(est)), rep("preds", nrow(sdat))))
  if (squared == TRUE){
    ddd$y <- ddd$y^2
    qul <- quu <- y@quantile^2
  }
  ymin <- min(ddd[,2][!is.na(ddd[,2])])
  ddd$x <- trans(ddd$x)
  pg <- function(x, y, group.number, ...){
    panel.abline(h=0, col="lightgrey")
    if (!is.null(qul)) panel.abline(h=qul, col="grey", lty=2)
    if (!is.null(quu)) panel.abline(h=quu, col="grey", lty=2)     
    if (group.number == 1){
      panel.segments(x, ymin, x, y, col="black", lwd=2)
    }
    if (group.number == 2){
      panel.segments(x, ymin, x, y, col="black", lwd=2, lty=3)
    }    
    if (group.number == 3){
      panel.lines(x, y, col="blue2", lwd=2)
    }    
  }
  xyplot(y ~ x | grp, group=typ, data=ddd, panel=panel.superpose, panel.groups=pg, scale=list(x=list(relation="free")), xlab=expression(beta), ylab="z", ...)
})




setMethod("plot", signature="mcptest", definition=function(x, y, trans=identity, alpha=NULL, order=TRUE, ...){
  est <- x@estimate
  pv <- x@pvalues
  if (order==TRUE){
    est <- est[order(pv, decreasing=TRUE)]
    pv <- pv[order(pv, decreasing=TRUE)]
  }
  color <- rep("lightgrey", length(est))
  if (!is.null(alpha)) color[pv < alpha] <- "orange"
  estr <- range(c(est, x@margin))
  estr[1] <- estr[1]-diff(estr)*0.05
  estr[2] <- estr[2]+diff(estr)*0.05
  par(mfrow=c(1,2))
  barplot(est, horiz=TRUE, las=2, xlab="estimates", xlim=estr, col=color, xpd=FALSE, ...)
  box()
  abline(v=0, lwd=2)  
  abline(v=x@margin, lwd=2, lty=2)
  barplot(rbind(pv), horiz=TRUE, xlim=c(0,1), xlab="p-values", col=color, ...)
  box()
  if (!is.null(alpha)) abline(v=alpha, lwd=2, lty=2)
})




setMethod("plot", signature=c(x="mcprofile", y="mcptest"), definition=function(x, y, squared=FALSE, trans=identity, ...){
  dp <- x@SRDP
  fs <- x@fsplines
  pspl <- lapply(1:length(fs), function(i){
    pseq <- seq(min(dp[[i]][,2]), max(dp[[i]][,2]), length=100)
    psp <- try(predict(fs[[i]], pseq), silent=TRUE)
    if (class(psp) == "try-error")  cbind(rep(NA, length(pseq)), rep(NA, length(pseq))) else cbind(psp$x, psp$y)
  })
  sx <- sapply(pspl, function(x) x[,1])
  sy <- sapply(pspl, function(x) x[,2])
  sdat <- data.frame(x=as.vector(sx), y=as.vector(sy))
  ffs <- as.factor(rep(names(dp), unlist(lapply(pspl, nrow))))

  est <- y@estimate
  zer <- rep(0, length(est))
  ffest <- as.factor(names(dp))
  margin <- y@margin

  tpspl <- lapply(1:length(fs), function(i){
    pseq <- seq(est[i], margin, length=100)
    psp <- try(predict(fs[[i]], pseq), silent=TRUE)
    if (class(psp) == "try-error")  cbind(rep(NA, length(pseq)), rep(NA, length(pseq))) else cbind(psp$x, psp$y)    
  })
  tsx <- sapply(tpspl, function(x) x[,1])
  tsy <- sapply(tpspl, function(x) x[,2])
  tsdat <- data.frame(x=as.vector(tsx), y=as.vector(tsy))
  tffs <- as.factor(rep(names(dp), unlist(lapply(tpspl, nrow))))
  
  ddd <- data.frame(x=c(sdat[,1], tsdat[,1], est), y=c(sdat[,2], tsdat[,2], zer), grp=c(as.character(ffs), as.character(tffs), as.character(ffest)), typ=c(rep("preds", nrow(sdat)), rep("tpreds", nrow(tsdat)), rep("ests", length(est))))
  if (squared == TRUE){
    ddd$y <- ddd$y^2
    qul <- quu <- y@quantile^2
  }
  ymin <- min(ddd[,2])
  ddd$x <- trans(ddd$x)
  pg <- function(x, y, group.number, ...){
    panel.abline(h=0, col="lightgrey")
    if (group.number == 2){
      panel.lines(x, y, col="blue2", lwd=2, lty=2)
    }
    if (group.number == 3){
      panel.lines(x, y, col="red2", lwd=3)
    }
    if (group.number == 1){
      panel.abline(v=x, col="black", lty=2)
    }
    panel.abline(v=margin, col="black", lty=2)
  }
  xyplot(y ~ x | grp, group=typ, data=ddd, panel=panel.superpose, panel.groups=pg, scale=list(x=list(relation="free")), xlab=expression(beta), ylab="z", ...)
})


setMethod("plot", signature=c(x="mcprofileRatio", y="mcptest"), definition=function(x, y, squared=FALSE, trans=identity, ...){
  dp <- x@SRDP
  fs <- x@fsplines
  pspl <- lapply(1:length(fs), function(i){
    pseq <- seq(min(dp[[i]][,2]), max(dp[[i]][,2]), length=100)
    psp <- try(predict(fs[[i]], pseq), silent=TRUE)
    if (class(psp) == "try-error")  cbind(rep(NA, length(pseq)), rep(NA, length(pseq))) else cbind(psp$x, psp$y)
  })
  sx <- sapply(pspl, function(x) x[,1])
  sy <- sapply(pspl, function(x) x[,2])
  sdat <- data.frame(x=as.vector(sx), y=as.vector(sy))
  ffs <- as.factor(rep(names(dp), unlist(lapply(pspl, nrow))))

  est <- y@estimate
  zer <- rep(0, length(est))
  ffest <- as.factor(names(dp))
  margin <- y@margin

  tpspl <- lapply(1:length(fs), function(i){
    pseq <- seq(est[i], margin, length=100)
    psp <- try(predict(fs[[i]], pseq), silent=TRUE)
    if (class(psp) == "try-error")  cbind(rep(NA, length(pseq)), rep(NA, length(pseq))) else cbind(psp$x, psp$y)    
  })
  tsx <- sapply(tpspl, function(x) x[,1])
  tsy <- sapply(tpspl, function(x) x[,2])
  tsdat <- data.frame(x=as.vector(tsx), y=as.vector(tsy))
  tffs <- as.factor(rep(names(dp), unlist(lapply(tpspl, nrow))))
  
  ddd <- data.frame(x=c(sdat[,1], tsdat[,1], est), y=c(sdat[,2], tsdat[,2], zer), grp=c(as.character(ffs), as.character(tffs), as.character(ffest)), typ=c(rep("preds", nrow(sdat)), rep("tpreds", nrow(tsdat)), rep("ests", length(est))))
  if (squared == TRUE){
    ddd$y <- ddd$y^2
    qul <- quu <- y@quantile^2
  }
  ymin <- min(ddd[,2])
  ddd$x <- trans(ddd$x)
  pg <- function(x, y, group.number, ...){
    panel.abline(h=0, col="lightgrey")
    if (group.number == 2){
      panel.lines(x, y, col="blue2", lwd=2, lty=2)
    }
    if (group.number == 3){
      panel.lines(x, y, col="red2", lwd=3)
    }
    if (group.number == 1){
      panel.abline(v=x, col="black", lty=2)
    }
    panel.abline(v=margin, col="black", lty=2)
  }
  xyplot(y ~ x | grp, group=typ, data=ddd, panel=panel.superpose, panel.groups=pg, scale=list(x=list(relation="free")), xlab=expression(beta), ylab="z", ...)
})




setMethod("plot", signature=c(x="mcprofile", y="mcprofile"), definition=function(x, y, squared=FALSE, trans=identity, ...){
  dp1 <- x@SRDP
  dp2 <- y@SRDP
  fs1 <- x@fsplines
  fs2 <- y@fsplines
  z1 <- sapply(dp1, function(x) x[,1])
  z2 <- sapply(dp2, function(x) x[,1])
  b1 <- sapply(dp1, function(x) x[,2])
  b2 <- sapply(dp2, function(x) x[,2])
  d1 <- data.frame(b=as.vector(unlist(b1)), z=as.vector(unlist(z1)))
  d2 <- data.frame(b=as.vector(unlist(b2)), z=as.vector(unlist(z2)))
  ff1 <- as.factor(rep(names(dp1), unlist(lapply(dp1, nrow))))
  ff2 <- as.factor(rep(names(dp2), unlist(lapply(dp2, nrow))))

  pspl1 <- lapply(1:length(fs1), function(i){
    pseq <- seq(min(dp1[[i]][,2]), max(dp1[[i]][,2]), length=100)
    psp <- try(predict(fs1[[i]], pseq), silent=TRUE)
    if (class(psp) == "try-error")  cbind(rep(NA, length(pseq)), rep(NA, length(pseq))) else cbind(psp$x, psp$y) 
  })
  sx1 <- sapply(pspl1, function(x) x[,1])
  sy1 <- sapply(pspl1, function(x) x[,2])
  sdat1 <- data.frame(x=as.vector(sx1), y=as.vector(sy1))
  ffs1 <- as.factor(rep(names(dp1), unlist(lapply(pspl1, nrow))))
   
  pspl2 <- lapply(1:length(fs2), function(i){
    pseq <- seq(min(dp2[[i]][,2]), max(dp2[[i]][,2]), length=100)
    psp <- try(predict(fs2[[i]], pseq), silent=TRUE)
    if (class(psp) == "try-error")  cbind(rep(NA, length(pseq)), rep(NA, length(pseq))) else cbind(psp$x, psp$y)
  })  
  sx2 <- sapply(pspl2, function(x) x[,1])
  sy2 <- sapply(pspl2, function(x) x[,2])
  sdat2 <- data.frame(x=as.vector(sx2), y=as.vector(sy2))
  ffs2 <- as.factor(rep(names(dp2), unlist(lapply(pspl2, nrow))))
  
  ddd <- data.frame(x=c(sdat1[,1],sdat2[,1]), y=c(sdat1[,2], sdat2[,2]), grp=c(as.character(ffs1), as.character(ffs2)), typ=c(rep("mcp1", nrow(sdat1)), rep("mcp2", nrow(sdat2))))
  if (squared == TRUE) ddd$y <- ddd$y^2
  ddd$x <- trans(ddd$x)
  pg <- function(x, y, group.number, ...){
    panel.abline(h=0, col="lightgrey")
    if (group.number == 1){
      panel.lines(x, y, col="blue2", lwd=2)
    }
    if (group.number == 2){
      panel.lines(x, y, col="red2", lwd=2, lty=2)
    }    
  }
  xyplot(y ~ x | grp, group=typ, data=ddd, panel=panel.superpose, panel.groups=pg, scale=list(x=list(relation="free")), xlab=expression(beta), ylab="z", key=list(space="top", text=list(c("profile 1", "profile 2")), lines=list(lty=c(1,2), col=c("blue2", "red2"), lwd=2)), ...)
})



