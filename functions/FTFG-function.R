###--------------------------------------------------------------###
###   Function for Functional Trend  Filtering on Graph (FTFG)   ###
###--------------------------------------------------------------###
## This code implements function FTFG


## package
library("fda.usc")


# get adjecency matrix
#adj <- hoge


# generate graph differnece operator matrix D from the adjecency matrix
#p <- 0; u <- NULL; v <- NULL; n <- nrow(adj)
#for (i in 1:(n-1)) {
#  for (j in (i + 1):n) {
#    if (adj[i, j] == 1) {
#      p <- p + 1
#      u <- c(u, i)
#      v <- c(v, j)
#    }
#  }
#}       
#m <- length(u); D <- matrix(0, m, n)
#for (p in 1:m) {
#  D[p, u[p]] <- 1
#  D[p, v[p]] <- -1
#}


# main function
FTFG <- function(fdata_, k, lam, L, ep1, ep2, D) {

### Carry out functional trend filtering 

## Arguments:
# fdata_ ...  functional data object of class "fdata"
# k      ...  the order of functional trend filtering (k=0,1,2)
# lam    ...  tuning parameter for reguralization
# L      ...  number of principal components
# ep1    ...  tolerance level
# ep2    ...  torelence level
# D      ...  graph differnece operator matrix

## Returns: A list with these named entries:
# beta_hat ... An estimated functional trend of class "fdata"
# KL       ... Karhunen–Loève estimation of trend of class "fdata"



  #  Check fdata_
  if (class(fdata_) != "fdata")  stop("First argument is not a functional data object of class fdata.")

  #  define a function named S_{\lambda} in the paper
  S_lam <- function(s,lam) max(0, (1-lam/norm(s,type="2")))*s

  #  preparations
  C  <- D  
  if (k>0) {
    for (i in 1:k) { 
      if (i%%2) C<-(t(D)%*%C) else C <-(D%*%C)
    }
  }
  n   <- nrow(fdata_)
  q   <- nrow(C)
  I   <- diag(n) 
  ea  <- diag(q)  
  eb  <- diag(L)
  rho <- 0.1

  #  get pc functions, mean function and their coefficients 
  mf   <- create.pc.basis(fdata_,l=1:L)$mean$data
  pcfs <- t(create.pc.basis(fdata_,l=1:L)$basis$data)
  X    <- create.pc.basis(fdata_,l=1:L)$x[,1:L]

  #  get optimal lambda for functional HP filter
  logspace <- function(start, end, n) exp(log(10)*seq(start, end, length.out=n))
  can <- logspace( -3, 3, 10 )
  MSE_optimal <- 10^10
  for (la in can) {
    fhpg <- FHPG(fdata(covid), k=1, la, L=5, D)$data
    MSE  <- 0
    for (j in 1:q) {
      MSE <- MSE+ norm((fhpg[,j]-covid[,j]), type="2")^2 /n/420
    }
    if (MSE < MSE_optimal) {
      MSE_optimal    <- MSE
      lambda_fhpg    <- la
    }  
  }

  #  set the initial values 
  j  <- 0
  b  <- FHPG(fdata(covid), k=1, lambda_fhpg, L=5, C)$coefs
  b0 <- b
  a  <- C%*%b    
  u  <- matrix(1, nrow=q, ncol=L)
  er <- 10^3

  #  itetration
  while (er > ep1) {

    #  update b
    for (l in 1:L) {
      b0[,l] <- b[,l]
      s1     <- colSums(u[,l]%*%C)
      s2     <- (eb[,l]%*%t(a))%*%C
      b[,l]  <- solve(I+rho*t(C)%*%C)%*%t(X[,l] - s1 + rho*s2)
    }

    #  update a 
    for (i in 1:q) {
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
    for (i in 1:q) {
      u[i,] <- u[i,]+rho*(ea[i,]%*%C%*%b-t(eb%*%a[i,]))
    }

    #  update er (sum of l2 error)
    er <- 0
    for (l in 1:L) {
      er <- er+norm((b[,l]-b0[,l]),type="2")/L
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
  return(result)
}


