glm_profiling <-
function(object, K, control, margin){
  est <- coefficients(object)
  Kest <- as.vector(K %*% est)
  sdKest <- sqrt(rbind(K) %*% vcov(object) %*% cbind(K))
  
  if (is.null(margin)){
    bup <- bdown <- numeric(length=control$maxsteps)
    bup[1] <- bdown[1] <- Kest
    bup[2] <- bup[1] + sdKest*0.001
    bdown[2] <- bdown[1] - sdKest*0.001
  } else {
    msequp <- seq(margin[1], margin[2], length=control$steps)
    mseqdown <- seq(margin[2], margin[1], length=control$steps)
    bup <- c(Kest, msequp[msequp >= Kest], msequp[control$steps])
    bdown <- c(Kest, mseqdown[mseqdown < Kest], mseqdown[control$steps])
  }

  zup <- numeric(length=length(bup))
  zdown <- numeric(length=length(bdown))

  pos <- abs(K) > 0
  if (sum(pos) > 1) nc <- sum(pos)-1 else nc <- 1
  parup <- matrix(nrow=length(bup), ncol=nc)
  pardown <- matrix(nrow=length(bdown), ncol=nc)
  if (sum(pos) > 1){
    w <- which(pos)[1]
    parup[1,] <- pardown[1,] <- est[-w][pos[-w]] 
  } else {
    parup[1,] <- pardown[1,] <- est[pos]
  }
  
  dispersion <- summary(object)$dispersion
  delta <- control$cutoff/control$steps

  start <- NULL
  for (i in 2:(length(bup)-1)){
    copt <- try(Kopt(object, K, bup[i], start), silent=TRUE)
    if (class(copt) == "try-error"){
      bup <- bup[1:(i-1)]
      zup <- zup[1:(i-1)]
      break
    }
    cdev <- copt$value
    start <- parup[i,] <- copt$par
    zup[i] <- sqrt(abs(cdev - deviance(object))/dispersion)
    if (is.null(margin)){
      if (zup[i]-zup[i-1] < .Machine$double.eps){
        bup[i+1] <- bup[i] + (bup[i]-bup[i-1])*5
      } else {
        bup[i+1] <- bup[i] + sign(bup[i]-bup[i-1])*(delta*((bup[i]-bup[i-1])/(zup[i]-zup[i-1])))
      }
      if (zup[i] > control$cutoff & i > 2){
        bup <- bup[1:i]
        zup <- zup[1:i]
        break
      }
      if (i == 2 & zup[i] > control$cutoff){
        bup[i+1] <- bup[1] + 0.001
      }
    }
  }
  bup <- bup[1:i]
  zup <- zup[1:i]
  parup <- parup[1:i,, drop=FALSE]
  
  start <- NULL
  for (i in 2:(length(bdown)-1)){
    copt <- try(Kopt(object, K, bdown[i], start), silent=TRUE)
    if (class(copt) == "try-error"){
      bdown <- bdown[1:(i-1)]
      zdown <- zdown[1:(i-1)]
      break
    }
    cdev <- copt$value
    start <- pardown[i,] <- copt$par
    zdown[i] <- -sqrt(abs(cdev - deviance(object))/dispersion)
    if (is.null(margin)){
      if ((zdown[i]-zdown[i-1]) > -1*.Machine$double.eps){
        bdown[i+1] <- bdown[i] + (bdown[i]-bdown[i-1])*5
      } else {
        bdown[i+1] <- bdown[i] + sign(bdown[i]-bdown[i-1])*(delta*((bdown[i]-bdown[i-1])/(zdown[i]-zdown[i-1])))
      }
      if (zdown[i] < -control$cutoff & i > 2){
        bdown <- bdown[1:i]
        zdown <- zdown[1:i]
        break
      }
      if (i == 2 & zdown[i] < -control$cutoff){
        bdown[i+1] <- bdown[1] - 0.001
      }
    }
  }
  bdown <- bdown[2:i]
  zdown <- zdown[2:i]
  pardown <- pardown[2:i,,drop=FALSE]
        
  b <- c(bdown, bup)
  z <- c(zdown, zup)
  optpar <- rbind(pardown, parup)[order(b),,drop=FALSE]
  list(stats=data.frame(b, z)[order(b),], param=optpar)
}

