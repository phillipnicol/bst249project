


#' @importFrom MCMCpack rdirichlet
#' @export
SpSlNLP <- function(y,X,phi,
                       warmup,iters,sigma2=1) {

  K <- length(phi)
  p <- ncol(X)
  n <- nrow(X)
  a <- rep(1,K+1)
  z <- rep(0, p)
  Pi <- rep(1,K+1)
  Pi[1] <- 0.9
  Pi[2:(K+1)] <- 0.1/K
  v0 <- 1
  beta <- as.vector(coef(glmnet(X,y, intercept =FALSE),s=0.1))[2:(p+1)]
  resid <- y - X%*%beta
  sigma2 <- sum(resid^2)/n
  sigma0 <- 1
  results <- list()
  results$z <- matrix(0,nrow=iters-warmup,ncol=p)
  results$pi <- matrix(0,nrow=iters-warmup,ncol=K+1)
  results$beta <- matrix(0,nrow=iters-warmup,ncol=p)

  #Hash table
  H <- new.env()

  cat("Starting sampler ... ... \n")

  Err <- diag(sigma2,nrow=n)
  for(i in 1:iters) {
    if(i < warmup) {
      cat("(Warmup) Iteration ", i, "/",iters,"\n")
    } else {
      cat("(Sampling) Iteration ", i, "/",iters,"\n")
    }

    for(j in 1:p) {
      log.prob <- rep(1,K+1)
      #Do p(z_j = 0)
      z[j] <- 0
      if(sum(z) == 0) {
        log.prob[1] <- log(Pi[1])+dmvnorm(y,sigma=Err,log=TRUE)
      } else {
        key <- get_key(z)
        val <- lookup(H,key)
        if(!is.null(val)) {
          log.prob[1] <- log(Pi[1])+val
        } else{
          Xz <- as.matrix(X[,z>0])
          z.cut <- z[z>0]
          phis <- phi[z.cut]
          y.null <- y[z==0]
          ML <- LA(Xz,y,phis,sigma2,p)$la
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
        } else {
          Xz <- as.matrix(X[,z>0])
          z.cut <- z[z>0]
          phis <- phi[z.cut]
          y.null <- y[z==0]
          ML <- LA(Xz,y,phis,sigma2,p)$la
          log.prob[k+1] <- log(Pi[k+1])+ML
          H[[key]] <- ML
        }
      }

      log.prob <- log.prob - max(log.prob)
      prob <- exp(log.prob)/sum(exp(log.prob))
      z[j] <- sample(1:(K+1),size=1,prob=prob)-1
    }

    if(i > warmup) {
      beta <- rep(0,p)
      Xz <- as.matrix(X[,z>0])
      if(ncol(Xz) > 0) {
        z.cut <- z[z>0]
        phis <- phi[z.cut]
        val <- LA(Xz,y,phis,sigma2,p)
        beta.cut <- rmvnorm(n=1,mean=val$beta.star,
                            sigma=solve(val$H))
        beta[z>0] <- beta.cut
      }
    }

    # Update counts
    for(j in 0:K) {
      a[j+1] <- 1+sum(z==j)
    }

    #Draw pi
    Pi <- rdirichlet(n=1, alpha=a)


    if(i > warmup) {
      #Save results
      results$z[i-warmup,] <- z
      results$beta[i-warmup,] <- beta
      results$pi[i-warmup,] <- Pi
    }
  }

  results
}

LA <- function(Xz,y,phis,sigma2,p) {
  pt <- ncol(Xz)
  la <- pt*sqrt(2)
  la <- la - (p/2)*log(2*pi*sigma2)
  la <- la + (pt/2)*log(2*pi)

  val <- LA.maximizer(Xz,y,phis,sigma2)
  beta.star <- val$beta
  H <- val$H

  f <- -1/(2*sigma2)*t(y-Xz%*%beta.star)%*%(y-Xz%*%beta.star)
  f <- f-sum(phis/(beta.star^2))-sum(beta.star^2/(2*phis))
  f <- as.vector(f)
  la <- la + 0.5*log(1/det(H)) + f
  val <- list()
  val$la <- la
  val$beta.star <- beta.star
  val$H <- H
  val
}

LA.maximizer <- function(Xz,y,phis,sigma2) {
  pt <- ncol(Xz)
  Lambda <- diag(1/(2*phis),nrow=pt)
  #initialize using ridge
  betat.start <- solve(t(Xz)%*%Xz+Lambda)%*%t(Xz)%*%y
  ixs <- which(abs(betat.start) < 0.01)
  betat.start[ixs] <- 0.01
  #betat.start <- rep(1,pt)
  betat.start <- as.vector(betat.start)
  pt <- ncol(Xz)
  tau <- 1
  beta <- betat.start
  while(tau > 10^{-4}) {
    #Compute gradient
    grad <- sigma2^{-1}*(t(Xz)%*%y - t(Xz)%*%Xz%*%beta)
    grad <- grad+2*phis/beta^3 - beta/phis

    #Compute HESSIAN
    H <- -sigma2^{-1}*(t(Xz)%*%Xz)-diag(6*phis/beta^4+1/phis,nrow=pt)
    new.diff <- as.vector(solve(H,-grad))
    beta.new <- new.diff+beta
    tau <- max(abs(beta-beta.new))
    beta <- beta.new
  }
  val <- list()
  val$beta.star <- beta
  val$H <- -H
  val
}







