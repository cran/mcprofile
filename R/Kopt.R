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
      if (any(!apply(object$family$linkinv(offs), 1, object$family$validmu) & apply(X[,neg, drop=FALSE], 1, function(x) all(x == 0)))){
        return(Inf)
      } else {
        cfm <- glm.fit(X[,neg, drop=FALSE],object$y,weights=object$prior.weights,start=est[neg],offset=offs,family=object$family,control=glm.control(maxit = 2000))
        return(cfm$deviance)
      }
    }
    w <- which(pos)[1]
    X <- model.matrix(object)
    if (is.null(start)) start <- est[-w][pos[-w]]
    return(optim(start, optd, method="BFGS"))
  } else {
    bopt <- numeric(length=length(K))
    bopt[(1:length(K))[pos]] <- b/K[pos]
    bopt[neg] <- est[neg]
    X <- model.matrix(object)
    offs <- X[,!neg,drop=FALSE] %*% bopt[!neg] + off
    cfm <- glm.fit(X[,neg, drop=FALSE],object$y,weights=object$prior.weights,start=est[neg],offset=offs,family=object$family,control=glm.control(maxit = 2000))
    return(list(par=b/K[pos], value=cfm$deviance))
  }
}

