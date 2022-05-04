
##p1 number of non null betas
##beta1 the absolute value of each
#' @export
simulateTest <- function(n,p,p1,beta1,phi.nlp,phi.normal) {
  #Number of methods x number of metrics
  results <- matrix(nrow=3,ncol=3)

  beta <- rep(0,p)
  beta[1:p1] <- sample(c(-1,1),size=p1,replace=TRUE)*beta1

  X <- matrix(rnorm(p*n),nrow=n,ncol=p)
  y <- X%*%beta+ rnorm(n,sd=1)
  y <- as.vector(y)
  #Center
  y <- y-mean(y)
  X <- scale(X,scale=TRUE)

  out.nlp <- SpSlNLP(y=y,X=X,phi=phi.nlp,warmup=50,iters=100)
  out.norm <- SpSlNormal(y=y,X=X,phi=phi.normal,warmup=50,iters=100)

  #Analyze out.nlp
  pm_beta <- colMeans(out.nlp$beta)
  rmse <- sqrt(mean((pm_beta - beta)^2))
  angle <- sum(pm_beta*beta)/(sqrt(sum(pm_beta^2)*sum(beta^2)))
  pip <- colMeans(out.nlp$z)

  results[1,1] <- rmse; results[1,2] <- angle
  results[1,3] <- mean(apply(out.nlp$z,1,function(x) mcc(preds=ifelse(x>0,1,0),actual=ifelse(beta>0,1,0))))

  #Analyze out.norm
  pm_beta <- colMeans(out.norm$beta)
  rmse <- sqrt(mean((pm_beta - beta)^2))
  angle <- sum(pm_beta*beta)/(sqrt(sum(pm_beta^2)*sum(beta^2)))
  pip <- colMeans(out.norm$z)

  results[2,1] <- rmse; results[2,2] <- angle
  results[2,3] <-  mean(apply(out.norm$z,1,function(x) mcc(preds=ifelse(x>0,1,0),actual=ifelse(beta>0,1,0))))


  return(results)
}






