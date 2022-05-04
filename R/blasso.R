

#' @importFrom statmod rinvgauss
#' @export

Blasso <- function(y,X,iters=2000,warmup=1000) {
  n <- nrow(X); p <- ncol(X)
  results <- list()
  results$beta <- matrix(0,nrow=iters-warmup,ncol=p)
  results$sigma2 <- rep(0,iters-warmup)
  results$lambda <- rep(0,iters-warmup)

  lambda = lambda_init(y,X)
  tau2 = rexp(p,rate = lambda^2)
  #beta = as.vector(lm(y~X-1)$coefficients)
  #beta = rep(0,p)
  beta = as.vector(coef(glmnet(X,y, intercept =FALSE),s=0.1))[2:(p+1)]
  sigma2 <- 1

  for(i in 1:iters) {
    if(i < warmup) {
      cat("(Warmup) Iteration ", i, "/",iters,"\n")
    } else {
      cat("(Sampling) Iteration ", i, "/",iters,"\n")
    }

    val <- bayes_lasso_emp_bayes(y,X,beta,sigma2,tau2,lambda)
    beta <- val$beta
    sigma2 <- val$sigma2
    lambda <- val$lambda
    tau2 <- val$tau2

    if(i > warmup) {
      results$beta[i-warmup,] <- beta
      results$sigma2[i-warmup] <- sigma2
      results$lambda[i-warmup] <- lambda
    }
  }

  return(results)
}











## Iteration of Bayesian LASSO Gibbs Sampler

### Assuming an empirical bayes updation step for lambda, the following two functions are for the Gibbs updation step, and for the initialization of lambda
bayes_lasso_emp_bayes <- function(y,X,beta,sigma2,tau2,lambda){
  n = dim(X)[1]
  p = dim(X)[2]
  D_inv = diag(1/tau2)
  A = t(X)%*%X+D_inv
  A_inv = solve(A)
  #A_inv = diag(tau2)
  mean = A_inv%*%t(X)%*%y
  var = sigma2*A_inv
  #bet_new = mvrnorm(1,mean,var)
  bet_new = t(rmvnorm(1,mean,var))
  err = y-X%*%bet_new
  scale_new = t(err)%*%err/2+t(bet_new)%*%D_inv%*%bet_new/2
  sig_new = invgamma::rinvgamma(1,(n-1)/2+(p/2),scale=1/scale_new)
  lamb_tau = lambda^2
  tau_new = rep(NA, length(tau2))
  for(i in 1:length(tau2)){
    mu_tau = sqrt(lamb_tau*sig_new/(bet_new[i]^2))
    a = rinvgauss(1,mean=mu_tau,shape=lamb_tau)
    tau_new[i] = 1/a
  }
  lamb_new = sqrt((2*p)/sum(tau_new))
  return(list(beta = bet_new,sigma2 = sig_new,tau2 = tau_new,lambda = lamb_new))
}



lambda_init <- function(y,X){
  n = dim(X)[1]
  p = dim(X)[2]
  beta = coef(glmnet(X,y),s=0.1)
  pred = predict(glmnet(X,y),X,s=0.1)
  sd = sd(y-pred)
  lambda = p*sd/sum(abs(beta))
  if(n>p){
    lm = lm(y~X)
    beta = lm$coefficients
    sd = sd(lm$residuals)
    lambda = p*sd/sum(abs(beta))
  }
  return(lambda)
}



### Instead, if we are using a Gamma hyperprior on lambda^2, with shape parameter r and scale parameter delta, we have the following Gibbs updation step
bayes_lasso_hyp_prior <- function(y,X,beta,sigma2,tau2,lambda,r,delta){
  n = dim(X)[1]
  p = dim(X)[2]
  D_inv = diag(1/tau2)
  A = t(X)%*%X+D_inv
  A_inv = solve(A)
  mean = A_inv%*%t(X)%*%y
  var = sigma2*A_inv
  bet_new = mvrnorm(1,mean,var)
  err = y-X%*%bet_new
  scale_new = t(err)%*%err/2+t(bet_new)%*%D_inv%*%bet_new/2
  #sig_new = rinvgamma(1,(n-1)/2+(p/2),scale=scale_new)
  sig_new = 1
  lamb_tau = lambda^2
  tau_new = rep(NA, length(tau2))
  for(i in 1:length(tau2)){
    mu_tau = sqrt(lamb_tau*sig_new/(bet_new[i]^2))
    a = rinvgauss(1,mean=mu_tau,shape=lamb_tau)
    tau_new[i] = 1/a
  }
  lamb_new = rgamma(1,p+r,sum(tau_new)/2+delta)
  return(bet_new,sig_new,tau_new,lamb_new)
}
