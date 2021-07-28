###--------------------------------------------------------###
###  Functions for selecting the tuning parameter of FDPC  ###
###--------------------------------------------------------###
## This code implements function FTFG


## package
library("freqdom.fda")


## function 
FDPC_select <- function(fdata_, q_candidates, L) {

 ### Carry out functional trend filtering 

 ## Arguments:
 # fdata_       ...  functional data object of class "fdata"
 # L            ...  number of principal components
 # q_candidates ...  candidates of tuning parameter q

 ## Returns:
 # q_optimal    ...  optimal lambda chosen from the candidates 


 #  Check fdata_
  if (class(fdata_) != "fdata")  stop("First argument is not a functional data object of class fdata.")


  d1 <- nrow(fdata_) 
  d2 <- floor(d1/5)
  error <- 10^20
  for (q in q_candidates) {
    error_new <- 0
    
    #  K fold cross-validation
    for (p in 1:d2) {
    
      #  test data
      fdata_test  <- fdata_$data[seq(p,d1,d2),]

      #  train data
      fdata_train <- fdata_$data[-seq(p,d1,d2),]
      
      fdpc <- t(fts.dpca(fd(t(fdata_train)),q=q, Ndpc = L)$Xhat$coefs)
      
      train_est <- fdpc

      for (i in 1:5) {
        if ((p+(d2-1)*i-d2>0) && (p+(d2-1)*i-(d2-1)<46)) {
          test_est <- (train_est[p+(d2-1)*i-d2,]+train_est[p+(d2-1)*i-(d2-1),])/2
          error_new <- error_new + norm((fdata_test[i,]-test_est),type="2")^2 
        }
      }
    }
    
    #  update optimal lambda and its error
    if (error >= error_new) {
      q_optimal <- q
      error <- error_new
    } else break
  } 
  return(q_optimal)
}
