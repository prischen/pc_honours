# Packages
library(hierNet)
library(glinternet)
library(lmSubsets) # exhaustive
library(ModelMetrics)
library(dplyr)

# Diabetes Data
data(diabetes, package = "lars")
x = diabetes$x
y = diabetes$y

# Repeated Cross Validation Error
set.seed(2020)

K = 10
m = 10
n = length(y)
reperror = rep(0,3*m)
reperror = matrix(reperror, nrow = 3)


t1 = Sys.time()

for (j in 1:m){
  set.seed(2020+j)
  error = rep(0,6*10)
  error = matrix(error, nrow = 6)
  for (r in 1:K){
    sets = sample(rep(1:K,n)[1:n],n)
    err = rep(0,K*6)
    err = matrix(err, nrow = 6)
    for (k in 1:K){
      # k = 7
      test_set = which(sets == k)
      train_set = (1:n)[-test_set]
      
      x_train = x[train_set,]
      y_train = y[train_set]
      x_test = x[test_set,]
      y_test = y[test_set]
      
      # hierNet
      fit_hier = hierNet.path(x_train, y_train, strong=TRUE, trace = 0) 
      fitcv_hier = hierNet.cv(fit,x_train,y_train,trace=0)
      fittedcv_hier =  hierNet.path(x_train, y_train, diagonal=FALSE, strong=TRUE, trace = 0, lamlist = fitcv$lamhat.1se)
      pred_hier = predict(fitted_cv, newx = x_test)
      err[1,k] = rmse(as.vector(pred), as.vector(y_test))
      
      # glinternet
      glin_fit = glinternet.cv(X = x_train, Y = y_train, numLevels = c(rep(1,dim(x_train)[2])))
      pred_glin = predict(glin_fit, X = x_test, lambdaType = c("lambdaHat1Std"))
      err[2,k] = rmse(as.vector(pred_glin), as.vector(y_test))
      
      
      # lmSubsets
      lm_model = lmIntSubsets(y_train~.*., data = data.frame(cbind(y_train,x_train)))
      aic = rep(NA,dim(lm_model$submodel)[1])
      for (i in 1:dim(lm_model$submodel)[1]){
        if (!is.na(lm_model$submodel$RSS[i])){
          mod = refit(lm_model, size = lm_model$submodel$SIZE[i], best = lm_model$submodel$BEST[i])
          aic[i] = AIC(mod)
        }
      }
      lm_model$aic = cbind(lm_model$submodel$SIZE, lm_model$submodel$BEST, aic)
      bestmod = refit(lm_model, size = lm_model$submodel$SIZE[which.min(aic)], best = lm_model$submodel$BEST[which.min(aic)])
      pred_lm = predict(bestmod, X = x_test)
      err[3,k] = rmse(as.vector(pred_lm), as.vector(y_test))
     
    }
    # vector with 6 entries: mean of errors for each package averaged over K iterations
    # placed inside a matrix: row = package, column = 
    error[,r] = apply(err, 1, mean)
    print(c("iteration:", r, "seed:", 2020 + j))
  }
  # matrix: 6 by 10
  # row: package, column: error for each seed
  reperror[,j] = apply(error, 1, mean)
}
t2 = Sys.time()
timetaken = t2-t1