## packages
library("ftsa")
library("freqdom.fda")
library("fda.usc")
source("FDPC-select.R")
source("FHP-function.R")
source("FTF-function.R")

n <- 50
L <- 5
rep <- 250


## define the true trend 

# scenerio2
true_trend2 <- matrix(0, nrow=n, ncol=120)
for (i in 1:n) {
   for (j in 1:120) {
      true_trend2[i,j] <- 40*sin((i+j)/12)
   }
}

# scenario3
bfs3 <- matrix(0, nrow=n, ncol=L)
for (l in 1:L) {
   bfs3[(10*(l-1)+1):(10*l), l] <- 1
}
b3 <- matrix(0, nrow=120, ncol=L)
for (i in 1:120) {
   b3[i,] <- c(40,0,40,0,40)  
}
true_trend3 <- bfs3%*%t(b3)

# scenario4
true_trend4 <- matrix(0, nrow=n, ncol=120)
for (i in 1:n) {
   for (j in 1:120) {
      true_trend4[i,j] <- j/3*(sin(4*i/n-2)+2*exp(-30*(4*i/n-2)^2))
   }
}


## set the noise variance
sigma <- 7


## preparation 
MSE1_tf  <- 0
MSE1_hp  <- 0
MSE1_pc  <- 0
MSE1_dpc <- 0
MSE2_tf  <- 0
MSE2_hp  <- 0
MSE2_pc  <- 0
MSE2_dpc <- 0
MSE3_tf  <- 0
MSE3_hp  <- 0
MSE3_pc  <- 0
MSE3_dpc <- 0
MSE4_tf  <- 0
MSE4_hp  <- 0
MSE4_pc  <- 0
MSE4_dpc <- 0




# the order of trend filtering
k <- 0


# the candidates of the tuning parameters
logspace     <- function(start, end, n) exp(log(10)*seq(start, end, length.out=n))
candidates   <- logspace( -3, 3, 60)
q.candidates <- 1:20


# iteration
for (m in 1:rep) {
  
   # generate noise  
  e <- matrix(0, nrow=n, ncol=120)
  for (i in 1:n) {
    e[i,] <- rnorm(120,mean=0,sd =sigma)
  } 
  
  ## scenario1
  x1 <-e
  
  if (m==1){
    lambda_optimal_ftf1 <- FTF_select(fdata(x1), k, L=5, ep1=0.01, ep2=0.0001, candidates)
    lambda_optimal_fhp1 <- FHP_select(fdata(x1), k, L=5, candidates)
    q_optimal1          <- FDPC_select(fdata(x1), q=q.candidates, L)
  }

  res1      <- FTF(fdata(x1), k, lam=lambda_optimal_ftf1, L=5, 0.0001,0.0001)
  res1_ftf  <- res1$ftf
 
  # Functional HP filter
  res1_fhp  <- FHP(fdata(x1), k, lam=lambda_optimal_fhp1, L=5)

  # The standard functional principle component method
  res1_KL  <- res1$KL

  # The dynamic functional principle component method
  res1_dpca <- fts.dpca(fd(t(x1)), q = q_optimal1, Ndpc=L)$Xhat


  for(i in 1:n){
    MSE1_tf  <- MSE1_tf + (norm(res1_ftf$data[i,] ,type="2"))^2 /n/120/rep 
    MSE1_hp  <- MSE1_hp + (norm(res1_fhp$data[i,] ,type="2"))^2 /n/120/rep 
    MSE1_pc  <- MSE1_pc + (norm(res1_KL$data[i,] ,type="2"))^2 /n/120/rep 
    MSE1_dpc <- MSE1_dpc + (norm(res1_dpca$coefs[,i],type="2"))^2 /n/120/rep 
  }



  ## scenario2
  x2 <- true_trend2 + e
  if (m==1){
    lambda_optimal_ftf2 <- FTF_select(fdata(x2), k, L=5, ep1=0.01, ep2=0.0001, candidates)
    lambda_optimal_fhp2 <- FHP_select(fdata(x2), k, L=5, candidates)
    q_optimal2          <- FDPC_select(fdata(x2), q=q.candidates, L)
  }

  res2      <- FTF(fdata(x2), k, lam=lambda_optimal_ftf2 , L=5, 0.0001,0.0001)
  res2_ftf  <- res2$ftf

  # Functional HP filter
  res2_fhp  <- FHP(fdata(x2), k, lam=lambda_optimal_fhp2, L=5)

  # The standard functional principle component method
  res2_KL  <- res2$KL

  # The dynamic functional principle component method
  res2_dpca <- fts.dpca(fd(t(x2)), q = q_optimal2, Ndpc=L)$Xhat

  for(i in 1:n){
    MSE2_tf  <- MSE2_tf + (norm(res2_ftf$data[i,]-true_trend2[i,] ,type="2"))^2 /n/120/rep 
    MSE2_hp  <- MSE2_hp + (norm(res2_fhp$data[i,]-true_trend2[i,] ,type="2"))^2 /n/120/rep 
    MSE2_pc  <- MSE2_pc + (norm(res2_KL$data[i,]-true_trend2[i,]  ,type="2"))^2 /n/120/rep 
    MSE2_dpc <- MSE2_dpc + (norm(res2_dpca$coefs[,i]-true_trend2[i,],type="2"))^2 /n/120/rep 
  }


  ## scenario3
  x3 <- true_trend3+e

  if (m==1){
    lambda_optimal_ftf3 <- FTF_select(fdata(x3), k=k, L=5, ep1=0.01, ep2=0.0001, candidates)
    lambda_optimal_fhp3 <- FHP_select(fdata(x3), k=k, L=5, candidates)
    q_optimal3          <- FDPC_select(fdata(x3), q=q.candidates, L)
  }

  # The proposed functional trend filtering
  res3      <- FTF(fdata(x3), k=k, lam=lambda_optimal_ftf3, L=5, 0.0001,0.0001)
  res3_ftf  <- res3$ftf

  # Functional HP filter
  res3_fhp  <- FHP(fdata(x3), k=k, lam=lambda_optimal_fhp3, L=5)

  # The standard functional principle component method
  res3_KL  <- res3$KL

  # The dynamic functional principle component method
  res3_dpca <- fts.dpca(fd(t(x3)), q = q_optimal3, Ndpc=L)$Xhat

  for(i in 1:n ){
    MSE3_tf  <- MSE3_tf + (norm(res3_ftf$data[i,] - true_trend3[i,],type="2"))^2 /n/120/rep 
    MSE3_hp  <- MSE3_hp + (norm(res3_fhp$data[i,] - true_trend3[i,],type="2"))^2 /n/120/rep 
    MSE3_pc  <- MSE3_pc + (norm(res3_KL$data[i,] - true_trend3[i,],type="2"))^2 /n/120/rep 
    MSE3_dpc <- MSE3_dpc + (norm(res3_dpca$coefs[,i] - true_trend3[i,],type="2"))^2 /n/120/rep   
  }



  ## scenario4
  x4 <- true_trend4 + e

  if (m==1){
    lambda_optimal_ftf4 <- FTF_select(fdata(x4), k, L=5, ep1=0.01, ep2=0.0001, candidates)
    lambda_optimal_fhp4 <- FHP_select(fdata(x4), k, L=5, candidates)
    q_optimal4          <- FDPC_select(fdata(x4), q=q.candidates, L)
  }


  # The proposed functional trend filtering
  res4      <- FTF(fdata(x4), k,lam=lambda_optimal_ftf4, L=5, 0.0001,0.0001)
  res4_ftf  <- res4$ftf

  # Functional HP filter
  res4_fhp  <- FHP(fdata(x4), k, lam=lambda_optimal_fhp4, L=5)

  # The standard functional principle component method
  res4_KL  <- res4$KL

  # The dynamic functional principle component method
  res4_dpca <- fts.dpca(fd(t(x4)), q = q_optimal4, Ndpc=L)$Xhat


  for(i in 1:n){
    MSE4_tf  <- MSE4_tf + (norm(res4_ftf$data[i,] - true_trend4[i,],type="2"))^2 /n/120/rep 
    MSE4_hp  <- MSE4_hp + (norm(res4_fhp$data[i,] - true_trend4[i,],type="2"))^2 /n/120/rep 
    MSE4_pc  <- MSE4_pc + (norm(res4_KL$data[i,] - true_trend4[i,],type="2"))^2 /n/120/rep 
    MSE4_dpc <- MSE4_dpc + (norm(res4_dpca$coefs[,i] - true_trend4[i,],type="2"))^2 /n/120/rep  
  }
}






