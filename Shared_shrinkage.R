#rm(list=ls())

library(corrplot)
library(RColorBrewer)
library(truncdist)
library(MASS)
library(TeachingDemos)
library(bayestestR)
library(dplyr)
library(ggplot2)
library(gtools)
library(truncnorm)
library(GIGrvg)

Shared_shrinkage <- function(Sim,MCMC,burn,THIN){
  
 
  data=data   ## Call the data
  X <- data$X  ## Exposure variable
  Y <- data$Y  ## Response 
  X_cov <- data$Cov  ## Covariate terms
  
  N = dim(X)[1]       ## Total number of individual
  p_cov <- dim(X_cov)[2] ## Number of covariate
  p_tot <- dim(X)[2]       ## Number of main effects
  n_interac <- nrow(combinations(p_tot, 2))  ## Number of interaction effects
  p <- p_tot+n_interac+p_cov+1;p    ## Total number of parameters in the model
  
  logX <- log(X+1)   ## log transformation & scaling
  sd_logx <- apply(logX,2,sd)
  logX_star <- apply(logX,2,function(x)(x/sd(x)))  
  # apply(logX_star,2,sd)
  
   
  rr <- combinations(p_tot, 2)        ## interaction combinations
  logZ <- array(NA,dim=c(N,n_interac))    ## design matrix for interaction terms
  for(j in 1:n_interac)
  {
    logZ[,j] <- logX_star[,rr[j,1]]*logX_star[,rr[j,2]]
  }
  sd_logz <-  apply(logZ,2,sd)
  logZ_star <- apply(logZ,2,function(x)(x/sd(x)))
  # apply(logZ_star,2,sd)
  
  X <- cbind(logX_star,logZ_star)   ## combining log-transformed main & Interaction effect
  X <- cbind(rep(1,N),X_cov,X)     ## adding covariates and intercept
  
  X <- as.matrix(X)
  
  (prev_rate <-  sum(Y)/N)
  v <- 7
  
  
  ## Initial value
  
  w <- rep(1,N)*0
  
  lambda <- rep(.01,N)
  
  beta <- rep(1,p)
  eta_true <- eta <- array(1,p_tot)
  theta_true <- theta <- array(1,n_interac)
  
  a <- 1
  b <- 1
  a1 <- a2 <- a3 <- a4 <- a5 <- 1
  b1 <- b2 <- b3 <- b4 <- b5 <- 1
  
  
  ## Beta ~ N(0,1/a*eta_j)   ## main
  ## gama ~ N(0,1/b*eta_j*eta_k* theta_jk)   ## interaction
  
  #MCMC <- 150000
  #burn <- 5000
  #THIN=10
  
  nthin <- MCMC/THIN
  BETA <- array(NA,dim=c(nthin,p))
  W <- array(NA,dim=c(nthin,N))
  ETA <- array(NA,dim=c(nthin,p_tot))
  THETA <- array(NA,dim=c(nthin,n_interac))
  LAMBDA <- array(NA,dim=c(nthin,N))
  A <- array(nthin)
  B <- array(nthin)
  
  store_count <- 0
  
  
  
  x <- X[,1:(p_tot+1+p_cov)]
  z <- X[,-(1:(p_tot+1+p_cov))]
  beta_main <- beta[1:(p_tot+1+p_cov)]
  beta_interac <- beta[-(1:(p_tot+1+p_cov))]
  start_time <- Sys.time()
  
  for (m in 1:MCMC){    ## MCMC chain
    
    for (i in 1:N){
      w[i] <- ifelse(Y[i]==1,rtruncnorm(1, a=0, b=Inf, mean = X[i,]%*%beta, sd = sqrt(1/lambda[i])),
                     rtruncnorm(1, a=-Inf, b=0, mean = X[i,]%*%beta, sd = sqrt(1/lambda[i])))
      
      lambda[i] <- rgamma(1,0.5*(v+1),0.5*(v+(w[i]-X[i,]%*%beta)^2))
      
    }
    
    
    eta_interac <- array(NA,n_interac)    ## design matrix for interaction terms
    for(i in 1:n_interac)
    {
      eta_interac[i] <- eta[rr[i,1]]*eta[rr[i,2]]
      
    }
    eta_mat = diag(c(rep(100,(p_cov+1)),1/(a*eta),1/(b*eta_interac*theta)))
    eta_mat_inv = diag(c(rep(1/100,(p_cov+1)),a*eta,b*eta_interac*theta))
    
    
    sig_inv <- diag(lambda)
    
    
    sig_main_inv <- diag(c(rep(0.01,(p_cov+1)),a*eta))
    sig_int_inv <- diag(b*eta_interac*theta)
    
    sig_main <- solve(t(x)%*%sig_inv%*%x+sig_main_inv)
    mu_main <- sig_main%*%(t(x)%*%sig_inv%*%(w-z%*%beta_interac))
    
    beta_main <- mvrnorm(1,mu_main,sig_main)
    
    
    sig_interac <- solve(t(z)%*%sig_inv%*%z+sig_int_inv)
    mu_interac <- sig_interac%*%(t(z)%*%sig_inv%*%(w-x%*%beta_main))
    
    beta_interac <- mvrnorm(1,mu_interac,sig_interac)
    beta <- c(beta_main,beta_interac)
    
    a <- rgamma(1,a3+0.5*p_tot,b3+0.5*sum(beta[(1+p_cov+1):(p_tot+1+p_cov)]^2*eta))
    
    for(i in 1:p_tot){
      eta_k <- eta[-i]
      temp  <- which(as.numeric(apply(rr,1,function(x)any(x==i)))==1)
      theta_temp <- theta[temp]
      beta_temp <- beta[temp+p_tot+1+p_cov]
      eta[i] <- rgamma(1,a1+0.5*p_tot,b1+0.5*(a*beta[i+1+p_cov]^2+b*sum(theta_temp*eta_k*beta_temp^2)))
    }
    
    
    temp2 <- array(NA,n_interac)
    for(i in 1:n_interac){
      temp2[i] <- eta[rr[i,1]]*eta[rr[i,2]]*theta[i]*beta[i+p_tot+1+p_cov]^2
    }
    
    b <- rgamma(1,a4+0.5*n_interac,b4+0.5*sum(temp2))
    
    for(i in 1:n_interac){theta[i] <- rgamma(1,a2+0.5,b2+0.5*b*(eta[rr[i,1]]*eta[rr[i,2]]*beta[i+p_tot+1+p_cov]^2))}
    
    ## storage
    if(m%%THIN==0){
      store_count <- store_count+1
      BETA[store_count,] <- beta
      ETA[store_count,] <- eta
      THETA[store_count,] <- theta
      A [store_count] <- a
      B [store_count] <- b
      LAMBDA[store_count,] <- lambda
      W[store_count,] <- w
    }
    # print(m)
  }
  
  end_time <- Sys.time()
  end_time - start_time
  
  out=list(BETA=BETA,ETA=ETA,THETA=THETA,A=A,B=B)
  return(out)
}


