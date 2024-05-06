###--------------------------------------------------###
###  Functions for Functional Trend Filtering (FTF)  ###
###--------------------------------------------------###
## This code implements two functions: FTF and FTF_select


## package
library("fda.usc")



# main function
FTF <- function(fdata_, k, lam, L, ep1, ep2) {

### Carry out functional trend filtering 

## Arguments:
# fdata_ ...  functional data object of class "fdata"
# k      ...  the order of functional trend filtering (k=0,1,2)
# lam    ...  tuning parameter for reguralization
# L      ...  number of principal components
# ep1    ...  tolerance level
# ep2    ...  torelence level

## Returns: A list with these named entries:
# beta_hat ... An estimated functional trend of class "fdata"
# KL       ... Karhunen–Loève estimation of trend of class "fdata"


 #  Check fdata_
  if (class(fdata_) != "fdata")  stop("First argument is not a functional data object of class fdata.")

  #  define a function named S_{\lambda} in the paper
  S_lam <- function(s,lam) max(0, (1-lam/norm(s,type="2")))*s    

  #  set up difference operator matrix D 
  n  <- nrow(fdata_)
  D <- matrix(0, n-1, n)
  for(i in 1:(n-1)){
    for(j in 1:n){
      if(i==j){ D[i, j] <- -1 }
      if((j-i)==1){ D[i, j] <- 1 }
    }
  }  
  C<-D
  if(k>0){
    for(i in 1:k){
      C <- D[1:(n-i-1), 1:(n-i)]%*%C
    }
  }

  
  #  preparations
  I   <- diag(n) 
  ea  <- diag(n-k-1)  
  eb  <- diag(L)  
  rho <- 0.1

  #  get pc functions, mean function and their coefficients 
  start <- proc.time()
  mf   <- create.pc.basis(fdata_,l=1:L)$mean$data 
  pcfs <- t(create.pc.basis(fdata_,l=1:L)$basis$data) 
  X    <- create.pc.basis(fdata_,l=1:L)$x[,1:L] 
  end <- proc.time()
  print("fPCA")
  print(summary(end - start))

  #  set the initial values 
  j  <- 0
  b  <- X
  b0 <- b
  a  <- C%*%b    
  u  <- matrix(1, nrow=n-k-1, ncol=L) 
  er <- 10^3
  iter_cnt <- 0

  #  itetration
  while (er > ep1) {
    iter_cnt <- iter_cnt + 1
    #  update b
    if (L==1){
      b0 <- b
      s1 <- colSums(t(u)%*%C)
      s2 <- (eb%*%t(a))%*%C
      b  <- solve(I+rho*t(C)%*%C)%*%t(X - s1 + rho*s2)
    } else {
      for (l in 1:L) {
        b0[,l] <- b[,l]
        s1     <- colSums(u[,l]%*%C)
        s2     <- (eb[,l]%*%t(a))%*%C
        b[,l]  <- solve(I+rho*t(C)%*%C)%*%t(X[,l] - s1 + rho*s2)
      }  
    }
    #  update a 
    end <- n-k-1
    for (i in 1:end) {
      s_new <- 1
      s     <- 1
      w_new <- a[i,]
      w_hat <- a[i,]
      w     <- a[i,]+10   
      while (norm((w-w_new),type="2") > ep2 ) {
        w<- w_new
        s<-s_new  
        w_new <- S_lam((1-rho)*w_hat+rho*(ea[i,]%*%C%*%b+u[i,]/rho),lam)
        s_new <- (1+sqrt(1+4*s**2))/2
        w_hat <- w_new+(s-1)*(w_new-w)/s_new
      }
      a[i,] <- c(w_new)
    }

    #  update u
    for (i in 1:end) {
      u[i,] <- u[i,]+rho*(ea[i,]%*%C%*%b-t(eb%*%a[i,]))
    }

    ## update er (sum of l2 error)
    er <- 0
    if (L==1){
      er <- er+norm((b-b0),type="2")
    } else {  
      for (l in 1:L) {
      er <- er+norm((b[,l]-b0[,l]),type="2")/L
      }
    }
  }

  KL <- create.pc.basis(fdata_, l=1:L)$mean
  KL$data <- t(t(X %*% t(pcfs)) + c(mf))
  KL$names$main <- "Karhunen–Loève"

  ftf <- create.pc.basis(fdata_, l=1:L)$mean
  ftf$data <- t(t(b %*% t(pcfs)) + c(mf))
  ftf$names$main <- "functional trend filtering"  

  result = list()
  result$ftf = ftf
  result$KL = KL
  print(iter_cnt)
  return(result)
}











FTF_select <- function(fdata_, k, L, ep1, ep2, candidates) {
### Select tuning parameters

## Arguments:
# fdata_     ...  functional data object of class "fdata"
# k          ...  the order of functional trend filtering (k=0,1,2)
# L          ...  number of principal components
# ep1        ...  tolerance level
# ep2        ...  torelence level
# candidates ...  candidates of tuning parameter lambda

## Returns:
# lambda_optimal ... optimal lambda chosen from the candidates 


  d1 <- nrow(fdata_) 
  d2 <- floor(d1/5)
  error <- 10^20

  for (lam in candidates) {
    error_new <- 0

    #  K fold cross-validation
    for (p in 1:(d2-1)) {
      #  test data
      fdata_test  <- fdata_$data[seq(p,d1,d2),]
      #  train data
      fdata_train <- fdata_$data[-seq(p,d1,d2),]
      
      #  apply functional trend filtering to the train data
      result <- FTF(fdata(fdata_train), k, lam, L, ep1, ep2)
      train_est <- t(result$ftf$data)

      for (i in 1:5) {
        if ((p+(d2-1)*i-d2>0) && (p+(d2-1)*i-(d2-1)<46)) {
          test_est <- (train_est[,p+(d2-1)*i-d2]+train_est[,p+(d2-1)*i-(d2-1)])/2
          error_new <- error_new + norm((fdata_test[i,]-test_est),type="2")^2 
        }
      }
    }
    #  update optimal lambda and its error
    if (error >= error_new) {
      lambda_optimal <- lam
      error <- error_new
    }  else break
  }
  return(lambda_optimal)
}




