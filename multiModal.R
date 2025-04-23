# A function to analyse multimodality (in mixtures of normal distributions)

library(data.table)
library(MASS)
library(mixtools)
library(tidyr)
library(ggplot2)
library(patchwork)

multiModal <- function(x) {
  fit1 <- fitdistr(x,"normal")
  fit2 <- normalmixEM(x,k=2,mu = range(x),lambda=c(.5,.5),sigma=c(1,1)*sd(x),maxit = 1e4)
  fit3 <- normalmixEM(x,k=3,mu = c(min(x),mean(x),max(x)),lambda = c(1,1,1)/3,sigma=c(1,1,1)*sd(x),maxit = 1e4)
  
  fit1dens <- function(x) {
    dnorm(x,mean=fit1$estimate[1],sd=fit1$estimate[2])
  }
  
  fit2dens <- function(x) {
      fit2$lambda[1] * dnorm(x,mean=fit2$mu[1],sd=fit2$sigma[1]) + 
      fit2$lambda[2] * dnorm(x,mean=fit2$mu[2],sd=fit2$sigma[2])
  }

  t2 <- optimize(fit2dens,interval = fit2$mu)$minimum
  
  fit3dens <- function(x) {
      fit3$lambda[1] * dnorm(x,mean=fit3$mu[1],sd=fit3$sigma[1]) + 
      fit3$lambda[2] * dnorm(x,mean=fit3$mu[2],sd=fit3$sigma[2]) +
      fit3$lambda[3] * dnorm(x,mean=fit3$mu[3],sd=fit3$sigma[3])
  }
  
  t3 <- c(optimize(fit3dens,interval = fit3$mu[1:2])$minimum,
          optimize(fit3dens,interval = fit3$mu[2:3])$minimum)
    
  param1 <- 2
  param2 <- 4 + 1
  param3 <- 6 + 2
  
  bic1 <- log(fit1$n) * param1 - 2*fit1$loglik
  bic2 <- log(length(fit2$x)) * param2 - 2*fit2$loglik
  bic3 <- log(length(fit3$x)) * param3 - 2*fit3$loglik
  
  prob2 <- 1-pchisq(fit2$loglik-fit1$loglik,param2-param1)
  prob3 <- 1-pchisq(fit3$loglik-fit2$loglik,param3-param2)
  
  best <- fit1
  if(prob2 < 0.01) best <- fit2
  if(prob3 < 0.01) best <- fit3
  
  message(sprintf("Model 2 vs Model 1 Likelihood Ratio Test: p=%1.2g",prob2))
  message(sprintf("Model 3 vs Model 2 Likelihood Ratio Test: p=%1.2g",prob3))
  
  d0 <- data.table(x=x,
                   group1=factor(1),
                   group2=factor(apply(fit2$posterior,1,which.max),levels=1:3),
                   group3=factor(apply(fit3$posterior,1,which.max),levels=1:3))

  n <- length(x)
  
  binwidth <- diff(range(x)) / 50
  
  xrange <- seq(from=min(x),to=max(x),l=1000)
 
  p1 <- ggplot(d0) + aes(x) +
    geom_histogram(aes(fill=group1),binwidth=binwidth) +
    geom_density(aes(y=after_stat(density * n * binwidth))) +
    geom_line(aes(x=x,y),data=data.table(x=xrange,y=fit1dens(xrange) * n * binwidth),lty=2) +
    annotate("label",x=min(x),y=Inf,label=sprintf("LL = %.1f\nBIC = %.1f",fit1$loglik,bic1),hjust=0,vjust=1) +
    ylab("Counts")
    
  p2 <- ggplot(d0) + aes(x) +
    geom_histogram(aes(fill=group2),binwidth=binwidth) +
    geom_density(aes(y=after_stat(density * n * binwidth))) +
    geom_line(aes(x=x,y),data=data.table(x=xrange,y=fit2dens(xrange) * n * binwidth),lty=2) +
    geom_vline(xintercept = t2,color="red") +
    annotate("label",x=min(x),y=Inf,label=sprintf("LL = %.1f\nBIC = %.1f\nM2 vs M1 p = %1.2g",fit2$loglik,bic2,prob2),hjust=0,vjust=1) +
    ylab("Counts")
  
  p3 <- ggplot(d0) + aes(x) +
    geom_histogram(aes(fill=group3),binwidth=binwidth) +
    geom_density(aes(y=after_stat(density * n * binwidth))) + 
    geom_line(aes(x=x,y),data=data.table(x=xrange,y=fit3dens(xrange) * n * binwidth),lty=2) +
    geom_vline(xintercept = t3,color="red") +
    annotate("label",x=min(x),y=Inf,label=sprintf("LL = %.1f\nBIC = %.1f\nM3 vs M2 p = %1.2g",fit3$loglik,bic3,prob3),hjust=0,vjust=1) +
    ylab("Counts")
  
    
  print(p1 + p2 + p3)
  
  invisible(list(fit1=fit1,fit2=fit2,fit3=fit3,f1d=fit1dens,f2d=fit2dens,f3d=fit3dens,t2=t2,t3=t3,best=best))
}

# run some examples

multiModal(rnorm(1000))
multiModal(c(rnorm(1000),rnorm(200,5,3)))
multiModal(c(rnorm(1000),rnorm(200,10,3),rnorm(300,102,2)))

           
