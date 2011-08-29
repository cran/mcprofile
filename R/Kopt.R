Kopt <-
function(object, K, b, start=NULL){
  off <- if (is.null(model.offset(object))) rep(0, length(object$y)) else model.offset(object)
  est <- coefficients(object)
  pos <- abs(K) > 0
  neg <- abs(K) == 0

  if (sum(pos) > 1){
    optd <- function(x){
      Ki <- -1*c(-b, K[-w][pos[-w]])/K[w]
      hv <- c(1,x)
      bo <- (Ki %*% hv)[1,1]
      bopt <- numeric(length=length(K))
      bopt[(1:length(K))[pos][-1]] <- x
      bopt[w] <- bo
      bopt[neg] <- est[neg]
      offs <- X[,!neg,drop=FALSE] %*% bopt[!neg] + off
      cfm <- glm.fit(X[,neg, drop=FALSE],object$y,weights=object$prior.weights,etastart=object$fitted,offset=offs,family=object$family,control=glm.control(maxit = 2000))
      cfm$deviance
    }
    w <- which(pos)[1]
    X <- model.matrix(object)
    if (sum(pos) == 2){
      sdp <- sqrt(diag(vcov(object))[pos][2])
      if (is.null(start)) start <- est[-w][pos[-w]]
      optres <- optimize(optd, interval=c(start-sdp*5, start+sdp*5))
      names(optres) <- c("par", "value")
      return(optres)
    } else {
      if (is.null(start)) start <- est[-w][pos[-w]]
      return(optim(start, optd, method="BFGS"))
    }
  } else {
    bopt <- numeric(length=length(K))
    bopt[(1:length(K))[pos]] <- b/K[pos]
    bopt[neg] <- est[neg]
    X <- model.matrix(object)
    offs <- X[,!neg,drop=FALSE] %*% bopt[!neg] + off
    cfm <- glm.fit(X[,neg, drop=FALSE],object$y,weights=object$prior.weights,etastart=object$fitted,offset=offs,family=object$family,control=glm.control(maxit = 2000))
    return(list(par=b, value=cfm$deviance))
  }
}

