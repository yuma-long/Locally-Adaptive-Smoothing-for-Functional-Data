###----------------------------------------------------------###
###   Function for Functional HP Filtering on Graph (FHPG)   ###
###----------------------------------------------------------###
## This code implements function FHPG


## package
library("fda.usc")


# get adjecency matrix
#adj <- hoge


# generate graph differnece operator matrix D from the adjecency matrix
#k <- 0; u <- NULL; v <- NULL; n <- nrow(adj)
#for (i in 1:(n-1)) {
#  for (j in (i + 1):n) {
#    if (adj[i, j] == 1) {
#      k <- k + 1
#      u <- c(u, i)
#      v <- c(v, j)
#    }
#  }
#}       
#m <- length(u); D <- matrix(0, m, n)
#for (k in 1:m) {
#  D[k, u[k]] <- 1
#  D[k, v[k]] <- -1
#}



# main function
FHPG <- function(fdata_, k, lam, L, D) {
### Carry out functional HP filtering on graph

## Arguments:
# fdata_ ...  functional data object of class "fdata"
# k      ...  the order of functional trend filtering (k=0,1,2)
# lam    ...  tuning parameter for reguralization
# L      ...  number of principal components
# D      ...  graph differnece operator matrix

## Returns:
# fhpg   ... functional HP estimation on graph of class "fdata"


  #  Check fdata_  
  if (class(fdata_) != "fdata")  stop("First argument is not a functional data object of class fdata.")
  
  #  preparation
  n  <- nrow(fdata_)
  C <- D
  if (k>0) {
    for (i in 1:k) { 
      if (i%%2) C<-(t(D)%*%C) else C <-(D%*%C)
    }
  }

  #  get pc functions, mean function and their coefficients 
  mf   <- create.pc.basis(fdata_,l=1:L)$mean$data 
  pcfs <- t(create.pc.basis(fdata_,l=1:L)$basis$data) 
  X    <- create.pc.basis(fdata_,l=1:L)$x[,1:L]
  I    <- diag(n) 
  
  fhpg            <- create.pc.basis(fdata_, l=1:L)$mean
  fhpg$coefs      <- solve(I+2*lam*t(C)%*%C)%*%X
  fhpg$data       <- t(pcfs %*% t(solve(I+2*lam*t(C)%*%C)%*%X)+c(mf))
  fhpg$names$main <- "functional HP filter on graph"
  return(fhpg)
}




