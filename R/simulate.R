
##p1 number of non null betas
##beta1 the absolute value of each
#' @importFrom horseshoe horseshoe
#' @importFrom mltools mcc
#' @importFrom mvtnorm rmvnorm
#' @export
simulateTest <- function(n,p,p1,beta1,phi.nlp,phi.normal) {
  #Number of methods x number of metrics
  results <- matrix(0,nrow=4,ncol=4)

  beta <- rep(0,p)
  beta[1:p1] <- sample(c(-1,1),size=p1,replace=TRUE)*beta1

  X <- matrix(rnorm(p*n),nrow=n,ncol=p)
  y <- X%*%beta+ rnorm(n,sd=1)
  print(dim(X))
  y <- as.vector(y)
  #Center
  y <- y-mean(y)
  X <- scale(X,scale=TRUE)

  out.nlp <- SpSlNLP(y=y,X=X,phi=phi.nlp,warmup=500,iters=1000)
  out.norm <- SpSlNormal(y=y,X=X,phi=phi.normal,warmup=500,iters=1000)
  out.bl <- Blasso(y=y,X=X,warmup=1000,iters=2000)
  out.horse = horseshoe::horseshoe(y,X,method.tau = "halfCauchy",method.sigma = "fixed",
                                   nmc=2000,burn=1000)

  #Analyze out.nlp
  pm_beta <- colMeans(out.nlp$beta)
  rmse <- sqrt(mean((pm_beta - beta)^2))
  angle <- sum(pm_beta*beta)/(sqrt(sum(pm_beta^2)*sum(beta^2)))
  pip <- colMeans(out.nlp$z)

  results[1,1] <- rmse; results[1,2] <- angle
  results[1,3] <- mean(apply(out.nlp$z,1,function(x) mcc(preds=ifelse(x>0,1,0),actual=ifelse(beta!=0,1,0))))
  results[1,4] <- mean(apply(out.nlp$z,1,function(x) {
    sum(x == 1 & beta != 0)/sum(beta!=0)
  }))#FN RATE


  #Analyze out.norm
  pm_beta <- colMeans(out.norm$beta)
  rmse <- sqrt(mean((pm_beta - beta)^2))
  angle <- sum(pm_beta*beta)/(sqrt(sum(pm_beta^2)*sum(beta^2)))
  pip <- colMeans(out.norm$z)

  results[2,1] <- rmse; results[2,2] <- angle
  results[2,3] <-  mean(apply(out.norm$z,1,function(x) mcc(preds=ifelse(x>0,1,0),actual=ifelse(beta!=0,1,0))))
  results[2,4] <- mean(apply(out.nlp$z,1,function(x) {
    sum(x == 1 & beta != 0)/sum(beta!=0)
  }))#FN RATE


  # Analyze out.bl
  pm_beta <- colMeans(out.bl$beta)
  rmse <- sqrt(mean((pm_beta - beta)^2))
  angle <- sum(pm_beta*beta)/(sqrt(sum(pm_beta^2)*sum(beta^2)))
  z <- apply(out.bl$beta,2,function(x) {
    lb <- quantile(x,0.025)
    ub <- quantile(x,0.975)
    if(0 < ub & 0 > lb) {
      0
    } else {
      1
    }
  })
  mcc <- mcc(preds=z,actual=ifelse(beta != 0,1,0))
  results[3,] <- c(rmse,angle,mcc,sum(z==1 & beta != 0)/sum(beta != 0))


  #Analyze horseshoe
  pm_beta <- out.horse$BetaHat
  betas <- out.horse$BetaSamples
  rmse <- sqrt(mean((pm_beta - beta)^2))
  angle <- sum(pm_beta*beta)/(sqrt(sum(pm_beta^2)*sum(beta^2)))
  z <- apply(betas,1,function(x) {
    lb <- quantile(x,0.025)
    ub <- quantile(x,0.975)
    if(0 < ub & 0 > lb) {
      0
    } else {
      1
    }
  })
  mcc <- mcc(preds=z,actual=ifelse(beta != 0,1,0))
  results[4,] <- c(rmse,angle,mcc,sum(z==1 & beta != 0)/sum(beta != 0))

  return(results)
}


## First simulation study with uncorrelated covariates

#' @export
FullSimulate1 <- function(reps,n,p,phi.nlp,phi.normal,beta1) {
  #Low dimensional setting
  #N and P

  #Number of non-nulls
  p1 <- c(10,20,30)
  nalg <- 4

  data <- matrix(0,nrow=0,ncol=5)
  for(i in 1:length(p1)) {
    mcc <- matrix(0,nrow=nalg,ncol=reps)
    fn <- matrix(0,nrow=nalg,ncol=reps)
    rmse <- matrix(0,nrow=nalg,ncol=reps)
    angle <- matrix(0,nrow=nalg,ncol=reps)
    for(j in 1:reps) {
      results <- simulateTest(n,p,p1[i],beta1,phi.nlp,phi.normal)
      rmse[,j] <- results[,1]
      angle[,j] <- results[,2]
      mcc[,j] <- results[,3]
      fn[,j] <- results[,4]
    }

    data <- rbind(data,c(rowMeans(rmse)[1],"SpSl-Normal",
                         p1[i],sd=sd(rmse[1,]),"RMSE"))
    data <- rbind(data,c(rowMeans(angle)[1],"SpSl-Normal",
                         p1[i],sd=sd(angle[1,]),"Cosine"))
    data <- rbind(data,c(rowMeans(mcc)[1],"SpSl-Normal",
                         p1[i],sd=sd(mcc[1,]),"MCC"))
    data <- rbind(data,c(rowMeans(fn)[1],"SpSl-Normal",
                         p1[i],sd=sd(fn[1,]),"FN"))

    data <- rbind(data,c(rowMeans(rmse)[2],"SpSl-NLP",
                         p1[i],sd=sd(rmse[2,]),"RMSE"))
    data <- rbind(data,c(rowMeans(angle)[2],"SpSl-NLP",
                         p1[i],sd=sd(angle[2,]),"Cosine"))
    data <- rbind(data,c(rowMeans(mcc)[2],"SpSl-NLP",
                         p1[i],sd=sd(mcc[2,]),"MCC"))
    data <- rbind(data,c(rowMeans(fn)[2],"SpSl-NLP",
                         p1[i],sd=sd(fn[2,]),"FN"))

    data <- rbind(data,c(rowMeans(rmse)[3],"BL",
                         p1[i],sd=sd(rmse[3,]),"RMSE"))
    data <- rbind(data,c(rowMeans(angle)[3],"BL",
                         p1[i],sd=sd(angle[3,]),"Cosine"))
    data <- rbind(data,c(rowMeans(mcc)[3],"BL",
                         p1[i],sd=sd(mcc[3,]),"MCC"))
    data <- rbind(data,c(rowMeans(fn)[3],"BL",
                         p1[i],sd=sd(fn[3,]),"FN"))

    data <- rbind(data,c(rowMeans(rmse)[4],"HS",
                         p1[i],sd=sd(rmse[4,]),"RMSE"))
    data <- rbind(data,c(rowMeans(angle)[4],"HS",
                         p1[i],sd=sd(angle[4,]),"Cosine"))
    data <- rbind(data,c(rowMeans(mcc)[4],"HS",
                         p1[i],sd=sd(mcc[4,]),"MCC"))
    data <- rbind(data,c(rowMeans(fn)[4],"HS",
                         p1[i],sd=sd(fn[4,]),"FN"))
  }

  data.LD <- data

  ## High dimensional setting Switch n and p
  n.new <- p
  p <- n
  n <- n.new

  data <- matrix(0,nrow=0,ncol=5)
  for(i in 1:length(p1)) {
    mcc <- matrix(0,nrow=nalg,ncol=reps)
    fn <- matrix(0,nrow=nalg,ncol=reps)
    rmse <- matrix(0,nrow=nalg,ncol=reps)
    angle <- matrix(0,nrow=nalg,ncol=reps)
    for(j in 1:reps) {
      results <- simulateTest(n,p,p1[i],beta1,phi.nlp,phi.normal)
      rmse[,j] <- results[,1]
      angle[,j] <- results[,2]
      mcc[,j] <- results[,3]
      fn[,j] <- results[,4]
    }

    data <- rbind(data,c(rowMeans(rmse)[1],"SpSl-Normal",
                         p1[i],sd=sd(rmse[1,]),"RMSE"))
    data <- rbind(data,c(rowMeans(angle)[1],"SpSl-Normal",
                         p1[i],sd=sd(angle[1,]),"Cosine"))
    data <- rbind(data,c(rowMeans(mcc)[1],"SpSl-Normal",
                         p1[i],sd=sd(mcc[1,]),"MCC"))
    data <- rbind(data,c(rowMeans(fn)[1],"SpSl-Normal",
                         p1[i],sd=sd(fn[1,]),"FN"))

    data <- rbind(data,c(rowMeans(rmse)[2],"SpSl-NLP",
                         p1[i],sd=sd(rmse[2,]),"RMSE"))
    data <- rbind(data,c(rowMeans(angle)[2],"SpSl-NLP",
                         p1[i],sd=sd(angle[2,]),"Cosine"))
    data <- rbind(data,c(rowMeans(mcc)[2],"SpSl-NLP",
                         p1[i],sd=sd(mcc[2,]),"MCC"))
    data <- rbind(data,c(rowMeans(fn)[2],"SpSl-NLP",
                         p1[i],sd=sd(fn[2,]),"FN"))

    data <- rbind(data,c(rowMeans(rmse)[3],"BL",
                         p1[i],sd=sd(rmse[3,]),"RMSE"))
    data <- rbind(data,c(rowMeans(angle)[3],"BL",
                         p1[i],sd=sd(angle[3,]),"Cosine"))
    data <- rbind(data,c(rowMeans(mcc)[3],"BL",
                         p1[i],sd=sd(mcc[3,]),"MCC"))
    data <- rbind(data,c(rowMeans(fn)[3],"BL",
                         p1[i],sd=sd(fn[3,]),"FN"))

    data <- rbind(data,c(rowMeans(rmse)[4],"HS",
                         p1[i],sd=sd(rmse[4,]),"RMSE"))
    data <- rbind(data,c(rowMeans(angle)[4],"HS",
                         p1[i],sd=sd(angle[4,]),"Cosine"))
    data <- rbind(data,c(rowMeans(mcc)[4],"HS",
                         p1[i],sd=sd(mcc[4,]),"MCC"))
    data <- rbind(data,c(rowMeans(fn)[4],"HS",
                         p1[i],sd=sd(fn[4,]),"FN"))
  }

  data.LD <- cbind(data.LD, "Low dimensional")
  data <- cbind(data, "High dimensional")
  data <- rbind(data.LD,data)

  return(data)
}



makePlot <- function(data) {
  data<-data.frame(metric=as.numeric(data[,1]),
                  Algorithm=data[,2],
                  S=as.numeric(data[,3]),
                  sd=as.numeric(data[,4]),
                  type=data[,5],
                  dimension=data[,6])

  p <- ggplot(data=data,aes(x=S,y=metric,color=Algorithm,
                            ymin=metric-sd,ymax=metric+sd))
  p <- p + geom_point()
  p <- p + geom_line()
  p <- p + geom_errorbar(width=0.1)
  p <- p + theme_linedraw()
  p <- p + xlab("# of Non-Nulls")
  p <- p + ylab("")
  p <- p + facet_grid(type ~ dimension)
}






