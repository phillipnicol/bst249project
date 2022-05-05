
#' @export
SpSlNormal <- function(y,X,phi,
                      warmup,iters) {

  K <- length(phi)
  n <- nrow(X)
  p <- ncol(X)
  a <- rep(1,K+1)
  z <- rep(0, p)
  v0 <- 1
  sigma0 <- 1
  sigma2 <- 1
  results <- list()
  results$z <- matrix(0,nrow=iters-warmup,ncol=p)
  results$pi <- matrix(0,nrow=iters-warmup,ncol=K+1)
  results$beta <- matrix(0,nrow=iters-warmup,ncol=p)
  results$sigma2 <- rep(0,iters-warmup)
  sigma2 <- sigma0
  Pi[1] <- 0.9
  Pi[2:(K+1)] <- 0.1/K
  #Hash table
  H <- new.env()

  cat("Starting sampler ... ... \n")


  for(i in 1:iters) {
    if(i < warmup) {
      cat("(Warmup) Iteration ", i, "/",iters,"\n")
    } else {
      cat("(Sampling) Iteration ", i, "/",iters,"\n")
    }

    for(j in 1:p) {
      log.prob <- rep(1,K+1)
      Err <- diag(sigma2,nrow=n)
      #Do p(z_j = 0)
      z[j] <- 0
      if(sum(z) == 0) {
        log.prob[1] <- log(Pi[1])+dmvnorm(y,sigma=Err,log=TRUE)
      } else {
        key <- get_key(z)
        val <- lookup(H,key)
        if(!is.null(val)) {
          log.prob[1] <- log(Pi[1]) + val
        } else {
          Xz <- as.matrix(X[,z>0])
          z.cut <- z[z>0]
          Psi.inv <- diag(1/phi[z.cut],nrow=length(z.cut))
          ML <- dmvnorm(y,sigma=Xz%*%Psi.inv%*%t(Xz)+Err,log=TRUE)
          log.prob[1] <- log(Pi[1])+ML
          H[[key]] <- ML
        }
      }

      for(k in 1:K) {
        z[j] <- k
        key <- get_key(z)
        val <- lookup(H,key)
        if(!is.null(val)) {
          log.prob[k+1] <- log(Pi[k+1])+val
        } else{
          Xz <- as.matrix(X[,z>0])
          z.cut <- z[z>0]
          Psi.inv <- diag(1/phi[z.cut],nrow=length(z.cut))
          ML <- dmvnorm(y,sigma=Xz%*%Psi.inv%*%t(Xz)+Err,log=TRUE)
          log.prob[k+1] <- log(Pi[k+1])+dmvnorm(y,sigma=Xz%*%Psi.inv%*%t(Xz)+Err,log=TRUE)
          H[[key]] <- ML
        }
      }

      log.prob <- log.prob - max(log.prob)
      prob <- exp(log.prob)/sum(exp(log.prob))
      z[j] <- sample(1:(K+1),size=1,prob=prob)-1
    }

    beta <- rep(0,p)
    Xz <- as.matrix(X[,z>0])
    if(ncol(Xz) > 0) {
      z.cut <- z[z>0]
      Psi.inv <- diag(1/phi[z.cut],nrow=length(z.cut))
      S <- Psi.inv+sigma2^{-1}*(t(Xz)%*%Xz)
      S.inv <- solve(S)
      beta.cut <- rmvnorm(n=1,mean=S.inv%*%(t(Xz)%*%y)/sigma2,
                            sigma=S.inv)
      beta <- rep(0,p)
      beta[z>0] <- beta.cut
    }

    # Update counts
    for(j in 0:K) {
      a[j+1] <- 1+sum(z==j)
    }
    Pi <- rdirichlet(n=1, alpha=a)

    ##Update sigma2
    fitted.mean <- X %*% beta
    SSR <- sum((y-fitted.mean)^2)
    shape <- 0.5*(v0+n)
    rate <- 0.5*(v0*sigma0+SSR)
    Gam <- rgamma(n=1,shape=shape,rate=rate)
    sigma2 <- 1/Gam


    if(i > warmup) {
      #Save results
      results$z[i-warmup,] <- z
      results$beta[i-warmup,] <- beta
      results$pi[i-warmup,] <- Pi
      results$sigma2[i-warmup] <- sigma2
    }
  }

  results
}

## Functions for hash table
lookup <- function(H, key) {
  val <- H[[key]]
  if(is.null(val)) {
    NULL
  } else {
    val
  }
}

get_key <- function(z) {
  paste(z, sep="", collapse="")
}
