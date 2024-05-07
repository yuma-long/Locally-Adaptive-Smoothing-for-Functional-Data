library(refund)
library(fda)
library(fds)
library(caret)

gas <- gasoline

wavelength <- as.numeric(gsub("nm", "", colnames(gas$NIR)))

plot(wavelength, gas$NIR[1,], type='l', ylim=range(gas$NIR), xlab="wave length", ylab="Value", main="Data")
for (i in 2:nrow(gas$NIR)) {
  lines(wavelength, gas$NIR[i,], col=rainbow(nrow(gas$NIR))[i])
}

# Split data into train and test

set.seed(42)
trainIndex <- createDataPartition(gasoline$octane, p = 0.8, list = FALSE)

trainSet <- gasoline[trainIndex,]
testSet <- gasoline[-trainIndex,]

model <- pffr(NIR ~ octane, yind = wavelength, data = list(NIR = trainSet$NIR, octane = trainSet$octane), 
              bs.yindex=list(bs="ps", k=50, m=c(2,1)))
prediction_result_linear <- predict(model, testSet)
class(prediction_result_linear)
rownames(prediction_result_linear) <- rownames(testSet)

# Create basis ------------------------------------------------------------

nbasis <-  50
order <-  4
basis <-  create.bspline.basis(rangeval=c(min(wavelength), max(wavelength)), nbasis = nbasis, norder = order)
plot(basis)
fd_obj <- Data2fd(argvals = wavelength, y=t(as.matrix(gas$NIR)), basisobj = basis)
plot(fd_obj)

new_matrix <- rbind(gas$octane, fd_obj$coefs) # first row: value of octane, second row~: coefficient of basis functions
rownames(new_matrix) <- c("octane", rownames(fd_obj$coefs))

trainval_matrix <- new_matrix[, trainIndex]
test_matrix <- new_matrix[, -trainIndex]
write.csv(trainval_matrix, "basis_coefs_trainval.csv", row.names=TRUE)
write.csv(test_matrix, "basis_coefs_test.csv", row.names = TRUE)


# Regression (Denoise) --------------------------------------------------------------

denoised_coefs <- read.csv("basis_coefs_denoised_test.csv", header=TRUE, row.names=1)
denoised_coefs <- data.frame(lapply(denoised_coefs, as.numeric))
denoised_coefs_matrix <- as.matrix(denoised_coefs) # denoised coefficient matrix without octane value

rownames(denoised_coefs_matrix) <- rownames(fd_obj$coefs)
colnames(denoised_coefs_matrix) <- colnames(test_matrix)

fd_obj_cp <- fd_obj
fd_obj_cp$coefs <- denoised_coefs_matrix
prediction_result_tv <- t(eval.fd(wavelength, fd_obj_cp))


# Result ------------------------------------------------------------------

linear_dif <- testSet$NIR - prediction_result_linear
tv_dif <- testSet$NIR - prediction_result_tv

sum_sq_linear_dif <- sum(linear_dif^2)
sum_sq_tv_dif <- sum(tv_dif^2)
