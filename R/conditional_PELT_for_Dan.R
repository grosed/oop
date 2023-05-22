
rm(list=ls())
library(changepoint.np)


#######################################################
########## Functions

### a sequence of quantiles
qnt.seq <- function(x, K=NULL, method="nonequally"){
  
  n <- length(x)
  if(is.null(K)){
    K <- round(4*log(n))
  }
  Q <- matrix(0, K, n+1)
  sorted.x = sort(x)
  yK = -1 + (2*(1:K)/K-1/K)
  c = -log(2*n-1)
  if(method=="equally"){
    pK = seq(0, 1, length.out = K+2)[-c(1, K+2)]
  }else{
    pK  = (1+exp(c*yK))^-1
  }
  
  qnt = rep(NA, K)
  for (i in 1:K){
    j  = as.integer((n-1)*pK[i] + 1)
    qnt[i] = sorted.x[j]
  }
  return(qnt)
  
}

### minus log-likelihood
minloglik <- function(x, Q, i, j){
  
  if (i>=j) {
    min.log.lik <- 0 
  }else{
    t <- j-i+1
    mi <- sum(x[i:j] < Q) + 0.5*sum(x[i:j]==Q)
    
    if(mi==0 | mi==t){
      min.log.lik <- 0 
    }else{
      min.log.lik <- -mi*log(mi/t)-(t-mi)*log(1-mi/t) 
    }
    
  }
  
  return(min.log.lik)  
  
}

### cost function for PELT
costMatrix = function(x, qnt) {
  
  N = length(x)
  K = length(qnt)
  C <- array(NA, dim=c(length(qnt), N, N))
  
  for(k in 1:K){
    for (i in 1:N){
      for (j in 1:N){
        C[k, i, j] = minloglik(x=x, Q=qnt[k], i=i, j=j)
      }
    }
  }
  
  return(C)
}

### cost function for PELT-conditional
costMatrix_c = function(x, qnt) {
  
  N = length(x)
  K = length(qnt)
  
  C = matrix(nr = N, nc = N)
  qnt0 = c(0, qnt, Inf)

  
  for (i in 1:N){
    
    for (j in 1:N){
      
      if (i>=j) {
        
        min.log.lik <- 0 
        
      }else{
        
        t <- j-i+1
        mi <- rep(NA, (length(qnt0)-1))
        f_mi <- rep(NA, (length(qnt0)-1))
        
        for(l in 1:(length(qnt0)-1)){
          
          mi[l] <- sum(((qnt0[l] < x[i:j])+(x[i:j] <= qnt0[l+1]))==2)
          
          if(mi[l]==0 | mi[l]==t){
            f_mi[l] <- 0 
          }else{
            f_mi[l] <- -mi[l]*log(mi[l]/t)
          }
          
        }
        
        min.log.lik <- sum(f_mi)
        
      }
      
      C[i, j] = min.log.lik
      
    }
  }
  
  return(C)
}

### Haynes et al (2017)
PELT = function (x, qnt, thr.c=3) {
  
  # initialization
  K = length(qnt)
  
  C = costMatrix(x, qnt)
  
  N = length(x)
  
  penalty = thr.c*log(N)
  
  Fvec = rep(0, N + 1)
  Fvec[1] = -penalty
  cngPoints = list(NULL) # initializing a null list
  Rvec = 0
  
  for (t in 1:N) {
    
    # compute all the costs up to time t
    meanC <- -2*log(2*N-1)*colMeans(-C[, (Rvec + 1), t, drop=F])
    partitionsCosts = Fvec[Rvec + 1] + meanC + penalty
    
    # get the new F(t) and its relative changepoint between R
    Fvec[(t) + 1] = min(partitionsCosts)
    cngPoint = Rvec[which.min(partitionsCosts)]
    cngPoints[[t+1]] = c(cngPoints[[cngPoint + 1]], cngPoint + 1)
    
    # make a vector of the same lenth to filter based on F(t)
    filter = (Fvec[Rvec + 1] + meanC) <= Fvec[(t) + 1]
    #print(c(filter))
    
    # append the new time to values in R that meet the condition
    Rvec = c(Rvec[filter], t)
    
  }
  
  return(cngPoints[[N]][-1])
  
}

### PELT-conditional (section 2.2)
PELT_c = function (x, qnt, thr.c=8) {
  
  # initialization
  K = length(qnt)
  
  C = costMatrix_c(x, qnt)
  
  N = length(x)
  
  penalty = K + 2*sqrt(K*thr.c*log(N)) + 2*thr.c*log(N)
  
  Fvec = rep(0, N + 1)
  Fvec[1] = -penalty
  cngPoints = list(NULL) # initializing a null list
  Rvec = 0
  
  for (t in 1:N) {
    
    # compute all the costs up to time t
    aggregatedC <- C[(Rvec + 1), t, drop=F]

    partitionsCosts = Fvec[Rvec + 1] + aggregatedC + penalty
    
    # get the new F(t) and its relative changepoint between R
    Fvec[(t) + 1] = min(partitionsCosts)
    cngPoint = Rvec[which.min(partitionsCosts)]
    cngPoints[[t+1]] = c(cngPoints[[cngPoint + 1]], cngPoint + 1)
    
    # make a vector of the same lenth to filter based on F(t)
    filter = (Fvec[Rvec + 1] + aggregatedC) <= Fvec[(t) + 1]
    
    # append the new time to values in R that meet the condition
    Rvec = c(Rvec[filter], t)
    
  }
  
  return(cngPoints[[N]][-1])
  
}

### PELT-max (section 2.3)
PELT_max = function (x, qnt, thr.c=3) {
  
  # initialization
  K = length(qnt)
  
  C = costMatrix(x, qnt)
  
  N = length(x)
  
  penalty = thr.c*log(N)
  #penalty = 1 + 2*sqrt(1*thr.c*log(N)) + 2*thr.c*log(N)
  
  Fvec = rep(0, N + 1)
  Fvec[1] = -penalty
  cngPoints = list(NULL) # initializing a null list
  Rvec = 0
  
  for (t in 1:N) {
    
    # compute all the costs up to time t
    maxC <- apply(2*C[, (Rvec + 1), t, drop=F], 2, max)
    partitionsCosts = Fvec[Rvec + 1] + maxC + penalty
    
    # get the new F(t) and its relative changepoint between R
    Fvec[(t) + 1] = min(partitionsCosts)
    cngPoint = Rvec[which.min(partitionsCosts)]
    cngPoints[[t+1]] = c(cngPoints[[cngPoint + 1]], cngPoint + 1)
    
    # make a vector of the same lenth to filter based on F(t)
    filter = (Fvec[Rvec + 1] + maxC) <= Fvec[(t) + 1]
    #print(c(filter))
    
    # append the new time to values in R that meet the condition
    Rvec = c(Rvec[filter], t)
    
  }
  
  return(cngPoints[[N]][-1])
  
}



#######################################################
########## Example



# single change point at (1/2)
int.length = 200
true.cpt <- c(int.length)+1
jump.size = 1
mu0 <- c(rep(0, 2*int.length))
mu <- c(rep(0, int.length), rep(jump.size, int.length))
n <- length(mu)
x = rnorm(n) + mu

plot(x)
lines(mu, type="s", col=7, lwd=2)

# quantiles
QC = c(2, 4, 8)
m = 2
K = min(round(n*0.9), round(QC[m]*log(n)))
qnt = qnt.seq(x, K=K)

# cp detection with three methods
PELT(x, qnt=qnt, thr.c=3)
PELT_c(x, qnt=qnt, thr.c=0.01)
PELT_max(x, qnt=qnt, thr.c=1)


