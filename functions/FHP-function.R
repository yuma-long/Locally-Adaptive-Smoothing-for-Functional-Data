###-----------------------------------------------###
###  Functions for Functional HP Filtering (FHP)  ###
###-----------------------------------------------###
## This code implements two functions: FHP and FHP_select


## package
library("fda.usc")



# main function
FHP <- function(fdata_, k, lam, L) {

### Carry out functional HP filtering 

## Arguments:
# fdata_ ...  functional data object of class "fdata"
# k      ...  the order of functional trend filtering (k=0,1,2)
# lam    ...  tuning parameter for reguralization
# L      ...  number of principal components

## Returns:
# fhp       ... functional HP estimation of class "fdata"


  #  Check fdata_
  if (class(fdata_) != "fdata")  stop("First argument is not a functional data object of class fdata.")
  
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

  #  get pc functions, mean function and their coefficients 
  mf   <- create.pc.basis(fdata_,l=1:L)$mean$data 
  pcfs <- t(create.pc.basis(fdata_,l=1:L)$basis$data) 
  X    <- create.pc.basis(fdata_,l=1:L)$x[,1:L] 
  I    <- diag(n) 

  fhp            <- create.pc.basis(fdata_, l=1:L)$mean
  fhp$coefs      <- solve(I+2*lam*t(C)%*%C)%*%X
  fhp$data       <- t(pcfs %*% t(solve(I+2*lam*t(C)%*%C)%*%X)+c(mf))
  fhp$names$main <- "functional HP filter"
  return(fhp)
}





FHP_select <- function(fdata_, k, L, candidates) {

### Select tuning parameters

## Arguments:
# fdata_     ...  functional data object of class "fdata"
# k          ...  the order of functional trend filtering (k=0,1,2)
# L          ...  number of principal components
# candidates ...  candidates of tuning parameter lambda

## Returns:
# lambda_optimal ... optimal lambda chosen from the candidates 

  d1 <- nrow(fdata_) 
  d2 <- floor(d1/5)
  error <- 10^20
  for (lam in candidates) {
    error_new <- 0
    #  K fold cross-validation
    for (p in 1:d2) {
      #  test data
      fdata_test  <- fdata_$data[seq(p,d1,d2),]

      #  train data
      fdata_train <- fdata_$data[-seq(p,d1,d2),]
      
      #  apply functional trend filtering to the train data
      result <- FHP(fdata(fdata_train), k, lam, L)
      train_est <- t(result$data)      

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
    } else break
  } 
  return(lambda_optimal)
}



