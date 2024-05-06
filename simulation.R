## packages
library("ftsa")
library("freqdom.fda")
library("fda.usc")
library(plotly)
source("./functions/FDPC-select.R")
source("./functions/FHP-function.R")
source("./functions/FTF-function.R")
source("./functions/FTFL1-function.R")

n <- 400 # number of time "t"
H <- 120 # number of data per func
L <- 5 # number of basis functions 
rep <- 250 # number of iterations


## define the true trend 

# scenerio2
true_trend2 <- matrix(0, nrow=n, ncol=H)
for (i in 1:n) {
   for (j in 1:H) {
      true_trend2[i,j] <- 40*sin((i+j)/12)
   }
}

# scenario3
bfs3 <- matrix(0, nrow=n, ncol=L)
for (l in 1:L) {
   bfs3[(10*(l-1)+1):(10*l), l] <- 1
}
b3 <- matrix(0, nrow=H, ncol=L)
for (i in 1:120) {
   b3[i,] <- c(40,0,40,0,40)  
}
true_trend3 <- bfs3%*%t(b3)

# scenario4
true_trend4 <- matrix(0, nrow=n, ncol=H)
for (i in 1:n) {
   for (j in 1:H) {
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
  e <- matrix(0, nrow=n, ncol=H)
  for (i in 1:n) {
    e[i,] <- rnorm(H,mean=0,sd =sigma)
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
    MSE1_tf  <- MSE1_tf + (norm(res1_ftf$data[i,] ,type="2"))^2 /n/H/rep 
    MSE1_hp  <- MSE1_hp + (norm(res1_fhp$data[i,] ,type="2"))^2 /n/H/rep 
    MSE1_pc  <- MSE1_pc + (norm(res1_KL$data[i,] ,type="2"))^2 /n/H/rep 
    MSE1_dpc <- MSE1_dpc + (norm(res1_dpca$coefs[,i],type="2"))^2 /n/H/rep 
  }



  ## scenario2
  x2 <- true_trend2 + e
  if (m==1){
    lambda_optimal_ftf2 <- FTF_select(fdata(x2), k, L=5, ep1=0.01, ep2=0.0001, candidates)
    lambda_optimal_fhp2 <- FHP_select(fdata(x2), k, L=5, candidates)
    q_optimal2          <- FDPC_select(fdata(x2), q=q.candidates, L)
  }
  
  start <- proc.time()
  res2      <- FTF(fdata(x2), k, lam=lambda_optimal_ftf2 , L=5, 0.0001,0.0001)
  end <- proc.time()
  print(summary(end - start))
  res2_ftf  <- res2$ftf

  # Functional HP filter
  res2_fhp  <- FHP(fdata(x2), k, lam=lambda_optimal_fhp2, L=5)

  # The standard functional principle component method
  res2_KL  <- res2$KL

  # The dynamic functional principle component method
  res2_dpca <- fts.dpca(fd(t(x2)), q = q_optimal2, Ndpc=L)$Xhat

  for(i in 1:n){
    MSE2_tf  <- MSE2_tf + (norm(res2_ftf$data[i,]-true_trend2[i,] ,type="2"))^2 /n/H/rep 
    MSE2_hp  <- MSE2_hp + (norm(res2_fhp$data[i,]-true_trend2[i,] ,type="2"))^2 /n/H/rep 
    MSE2_pc  <- MSE2_pc + (norm(res2_KL$data[i,]-true_trend2[i,]  ,type="2"))^2 /n/H/rep 
    MSE2_dpc <- MSE2_dpc + (norm(res2_dpca$coefs[,i]-true_trend2[i,],type="2"))^2 /n/H/rep 
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
    MSE3_tf  <- MSE3_tf + (norm(res3_ftf$data[i,] - true_trend3[i,],type="2"))^2 /n/H/rep 
    MSE3_hp  <- MSE3_hp + (norm(res3_fhp$data[i,] - true_trend3[i,],type="2"))^2 /n/H/rep 
    MSE3_pc  <- MSE3_pc + (norm(res3_KL$data[i,] - true_trend3[i,],type="2"))^2 /n/H/rep 
    MSE3_dpc <- MSE3_dpc + (norm(res3_dpca$coefs[,i] - true_trend3[i,],type="2"))^2 /n/H/rep   
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
    MSE4_tf  <- MSE4_tf + (norm(res4_ftf$data[i,] - true_trend4[i,],type="2"))^2 /n/H/rep 
    MSE4_hp  <- MSE4_hp + (norm(res4_fhp$data[i,] - true_trend4[i,],type="2"))^2 /n/H/rep 
    MSE4_pc  <- MSE4_pc + (norm(res4_KL$data[i,] - true_trend4[i,],type="2"))^2 /n/H/rep 
    MSE4_dpc <- MSE4_dpc + (norm(res4_dpca$coefs[,i] - true_trend4[i,],type="2"))^2 /n/H/rep  
  }
}

res1_ftf_data <- res1_ftf&data
res2_ftf_data <- res2_ftf$data
res3_ftf_data <- res3_ftf$data
res4_ftf_data <- res4_ftf$data

# fig <- plot_ly(x = rep(1:50, each=120), y = rep(1:120, times=50), z = c(t(res4_ftf_data)),
#                marker = list(size = 2, color = 'rgba(132, 126, 243, .7)'))
fig <- plot_ly(type='scatter3d', mode='markers')
fig <- fig %>% layout(scene = list(xaxis = list(title = "t"),
                                   yaxis = list(title = "X"),
                                   zaxis = list(title = "value"),
                                   aspectratio = list(x = 1, y = 1, z = 1))
                      )
fig <- fig %>% add_trace(x = rep(1:n, each=H), y = rep(1:H, times=n), z = c(t(true_trend4+e)),
                         marker = list(size = 1, color = '#000000'), type='scatter3d', mode='markers')
fig <- fig %>% add_trace(x = rep(1:n, each=H), y = rep(1:H, times=n), z = c(t(res4_ftf_data)),
                         marker = list(size = 2, color = 'rgba(240, 60, 60, .6)'), type='scatter3d', mode='markers')
fig <- fig %>% add_trace(x = rep(1:n, each=H), y = rep(1:H, times=n), z = c(t(res4_ftfl1_data)),
                         marker = list(size = 2, color = 'rgba(60, 60, 240, .6)'), type='scatter3d', mode='markers')
fig



# ftf_l1 ------------------------------------------------------------------

lambda_optimal_ftfl1_2 <- FTFL1_select(fdata(x2), k, L=5, ep1=0.01, ep2=0.0001, candidates)
start <- proc.time()
res2_l1 <- FTFL1(fdata(x2), k, lam = lambda_optimal_ftfl1_2, L = 5, ep1 = 0.0001, ep2 = 0.0001)
end <- proc.time()
print(summary(end - start))
res2_ftfl1 <- res2_l1$ftf
res2_ftfl1_data <- res2_ftfl1$data

lambda_optimal_ftfl1_3 <- FTFL1_select(fdata(x3), k, L=5, ep1=0.01, ep2=0.0001, candidates)
res3_l1 <- FTFL1(fdata(x3), k, lam = lambda_optimal_ftfl1_2, L = 5, ep1 = 0.0001, ep2 = 0.0001)
res3_ftfl1 <- res3_l1$ftf
res3_ftfl1_data <- res3_ftfl1$data

lambda_optimal_ftfl1_4 <- FTFL1_select(fdata(x4), k, L=5, ep1=0.01, ep2=0.0001, candidates)
res4_l1 <- FTFL1(fdata(x4), k, lam = lambda_optimal_ftfl1_2, L = 5, ep1 = 0.0001, ep2 = 0.0001)
res4_ftfl1 <- res4_l1$ftf
res4_ftfl1_data <- res4_ftfl1$data


# plot 2d -----------------------------------------------------------------

x_sample <-40

true_matrix <- true_trend4
raw_data_matrix <- true_trend4+e
ftf_matrix <- res4_ftf_data
ftfl1_matrix <- res4_ftfl1_data

true_vector_t <- true_matrix[,x_sample]
raw_data_vector_t <- raw_data_matrix[, x_sample]
ftf_vector_t <- ftf_matrix[,x_sample]
ftfl1_vector_t <- ftfl1_matrix[,x_sample]

x_range = 1:length(raw_data_vector_t)

fig1 <- plot_ly(type='scatter')
fig1 <- fig1 %>% add_trace(x = x_range, y = true_vector_t,
                         line = list(color = '#ABAB24', dash = "dash"), mode = "lines", name = "True")
fig1 <- fig1 %>% add_trace(x = x_range, y = raw_data_vector_t, 
                         marker = list(size = 5, color = '#000000'), mode = "markers", name = "Data")
fig1 <- fig1 %>% add_trace(x = x_range, y = ftf_vector_t,
                         marker = list(size = 5, color = 'rgba(240, 60, 60, .6)'), mode = "lines+markers", name = "FTF")
fig1 <- fig1 %>% add_trace(x = x_range, y = ftfl1_vector_t,
                         marker = list(size = 5, color = 'rgba(30, 30, 240, .6)'), mode = "lines+markers", name = "FTFL1")
fig1

t_sample = 10

true_vector_x <- true_matrix[t_sample,]
raw_data_vector_x <- raw_data_matrix[t_sample,]
ftf_vector_x <- ftf_matrix[t_sample,]
ftfl1_vector_x <- ftfl1_matrix[t_sample,]

t_range = 1:length(raw_data_vector_x)

fig2 <- plot_ly(type='scatter')
fig2 <- fig2 %>% add_trace(x = t_range, y = true_vector_x,
                         line = list(color = '#ABAB24', dash = "dash"), mode = "lines", name = "True")
fig2 <- fig2 %>% add_trace(x = t_range, y = raw_data_vector_x, 
                         marker = list(size = 5, color = '#000000'), mode = "markers", name = "Data")
fig2 <- fig2 %>% add_trace(x = t_range, y = ftf_vector_x,
                         marker = list(size = 5, color = 'rgba(240, 60, 60, .6)'), mode = "lines+markers", name = "FTF")
fig2 <- fig2 %>% add_trace(x = t_range, y = ftfl1_vector_x,
                         marker = list(size = 5, color = 'rgba(30, 30, 240, .6)'), mode = "lines+markers", name = "FTFL1")
fig2


# principal component functions -------------------------------------------

pcfs <- res2$pcfs
pcfs_df <- data.frame(pcfs)

fig3 <- plot_ly(pcfs_df, type = 'scatter', mode = 'lines')
fig3 <- fig3 %>% add_trace(y = ~PC1, name = 'PC1', mode = "lines")
fig3 <- fig3 %>% add_trace(y = ~PC2, name = 'PC2', mode = "lines")
fig3 <- fig3 %>% add_trace(y = ~PC3, name = 'PC3', mode = "lines")
fig3 <- fig3 %>% add_trace(y = ~PC4, name = 'PC4', mode = "lines")
fig3 <- fig3 %>% add_trace(y = ~PC5, name = 'PC5', mode = "lines")
fig3
