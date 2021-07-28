###---------------------------------------------------------------------------------------------###
###   estimate trend of transition of COVID-19 cases per one million population by prefecture   ###
###---------------------------------------------------------------------------------------------###

library("fda.usc")
library("NipponMap")
library("RColorBrewer")
source("FHPG-function.R")
source("FTFG-function.R")


# get the adjecency matrix of prefectures in Japan
adj <- read.csv(file="AdjacencyMatrix_JPN.csv")


# get the number of infected people per million in each prefecture from 1/16/2020 to 3/9/2021
covid <- as.matrix(read.csv("CovidPerPopulation.csv"))


# generate graph differnece operator matrix D from adjecency matrix
p <- 0; u <- NULL; v <- NULL; n <- nrow(adj)
for (i in 1:(n-1))   for (j in (i + 1):n)  if (adj[i, j] == 1) {
  p <- p + 1; u <- c(u, i); v <- c(v, j)
}   
m <- length(u); D <- matrix(0, m, n)
for (p in 1:m) { D[p, u[p]] <- 1; D[p, v[p]] <- -1}



###-------------------------------###
###  Apply FHPG to COVID-19 data  ###
###-------------------------------###

# decide the order of the methods
k <- 2

# prepare candidates of lambda
logspace <- function(start, end, n) exp(log(10)*seq(start, end, length.out=n))
candidates <- logspace( -3, 3, 60 )

# select lambda, which minimizes MSE 
MSE_optimal <- 10^10
for (lam in candidates) {
  fhpg <- FHPG(fdata(covid), k=k, lam=lam, L=5, D)$data
  MSE<-0
  for (j in 1:n) {
      MSE <- MSE+ norm((fhpg[,j]-covid[,j]), type="2")^2 /n/420
  }
  if (MSE < MSE_optimal) {
    MSE_optimal    <- MSE
    lambda_optimal <- lam
  }  
}

# JapanPrefMap
fhpg <- FHPG(fdata(covid), k=k, lam=lambda_optimal, L=5, D)$data[,395]
breaks <- c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 1000)
fhpg_cols <- as.character(cut(fhpg, breaks=breaks, labels=brewer.pal(9, "YlOrRd"), include.lowest=F))
JapanPrefMap(col = fhpg_cols, main = "FHP (k=1)") 



###-------------------------------###
###  Apply FTFG to COVID-19 data  ###
###-------------------------------###

# prepare candidates of lambda
logspace <- function(start, end, n) exp(log(10)*seq(start, end, length.out=n))
candidates <- logspace( -3, 3, 60 )

# select lambda, which minimizes MSE 
MSE_optimal <- 10^10
for (lam in candidates) {
  ftfg <- FTFG(fdata(covid), k=k, lam=lam, L=5, 0.001, 10^(-3), D)$ftf$data
  MSE  <- 0
  for (j in 1:n) {
      MSE <- MSE+ norm((ftfg[,j]-covid[,j]), type="2")^2 /n/420
  }
  if (MSE < MSE_optimal) {
    MSE_optimal    <- MSE
    lambda_optimal <- lam
  }  
}

# JapanPrefMap
ftfg <- FTFG(fdata(covid), k=k, lam=lambda_optimal, L=5, 10^(-3), 10^(-6), D)$ftf$data[,395]
breaks <- c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 1000)
ftfg_cols <- as.character(cut(round(ftfg), breaks=breaks, labels=brewer.pal(9, "YlOrRd"), include.lowest=F))
JapanPrefMap(col = ftfg_cols, main = "FTF (k=1)") 



