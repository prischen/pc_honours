# Packages
library(hierNet)
library(glinternet)
library(lmSubsets) # exhaustive
library(ModelMetrics)
library(dplyr)
library(lars)
library(doParallel)
library(foreach)
library(tidyr)
library(mvtnorm)

# Save Path
path = paste0(getwd(),"/HierProcessing/Results/",Sys.Date(),"/")
dir.create(path)

# Parameters
p = 10
n = 10*p

# Setting up the Data
Sigma = diag(1,p) 
rho = 0.2
for(i in 1:p){
  for(j in 1:p){
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = rmvnorm(n = n, mean = rep(0,p), sigma = Sigma) # Generating the regressors
colnames(X) = paste0("X",1:p)
X = data.frame(X)

# Data
data = model.matrix(~ 0 + .*., data.frame(X))

# Main Effects Only
y1 = 5*X$X1 + 3*X$X3 -2*X$X6 + X$X10
data1 = data.frame(y = y1, data)
# Main Effects + 2 Strong Hierarchy Interactions
y2 = 5*X$X1 + 3*X$X3 -2*X$X6 + X$X10 + 2*X$X6*X$X1 - 4*X$X3*X$X10
data2 = data.frame(y = y2, data)
# Main Effects + 2 Non - ME Interactions
y3 = 5*X$X1 + 3*X$X3 -2*X$X6 + X$X10 + 2*X$X2*X$X4 - 4*X$X7*X$X8
data3 = data.frame(y = y3, data)
# Main Effects + 2 ME - Non - ME Interactions
y4 = 5*X$X1 + 3*X$X3 -2*X$X6 + X$X10 + 2*X$X2*X$X1 - 4*X$X7*X$X10
data4 = data.frame(y = y4, data)


foreach(j=1:4) %dopar% {
  
  models = rep(0,3)
  
  y = get(paste0("y",j))  
  # hierNet
  
  s_time = Sys.time()
  fit_hier = hierNet.path(X, y, strong=TRUE, trace = 0) 
  fitcv_hier = hierNet.cv(fit_hier,X,y,trace=0)
  fittedcv_hier =  hierNet.path(X, y, diagonal=FALSE, strong=TRUE, trace = 0, lamlist = fitcv_hier$lamhat.1se)
  
  hier_features = names(fittedcv_hier$mx)
  hier_coeff = fittedcv_hier$bp - fittedcv_hier$bn
  intmat = fittedcv_hier$th
  hier_int = c()
  for(g in 1:(length(hier_features))){
    for(h in g:length(hier_features)){
      if (abs(intmat[g,h,1]) >= tol){
        intt = paste0(hier_features[g],":",hier_features[h])
        hier_int = c(hier_int, intt)
      }
    }
  }
  
  models[1] = paste0(c(hier_features[abs(hier_coeff) >= tol],hier_int), collapse = ", ")
  
  # glinternet
  
  glin_fit = glinternet.cv(X = X, Y = y, numLevels = c(rep(1,dim(X)[2])))
  glin_mod = 
    glin_features = names(X)
  me_cont = glin_mod$activeSet[[2]]$cont
  glin_me = glin_features[me_cont]
  int_cont = glin_mod$activeSet[[2]]$contcont
  glin_int = c()
  
  if (length(int_cont) > 0){
    for(q in 1:dim(int_cont)[1]){
      intt = paste0(glin_features[int_cont[q,1]],":",glin_features[int_cont[q,2]])
      glin_int = c(glin_int, intt)
    }
  }
  
  glin_mod$activeSet[[3]] = convert_glinternet_features(glin_mod$activeSet, cont_features = glin_features)
  
  models[2] = paste0(c(glin_me,glin_int), collapse = ", ")
  
  # lmSubsets w/ BIC
  lm_mod = lmIntSubsets(y~.*., data = data.frame(cbind(y,X)), nbest = 1000, nmax = 30)
  bestmod = lmSelect(lm_mod, penalty = "BIC")
  mod = refit(bestmod, data = data.frame(cbind(y,X)))
  models[3] = paste0(variable.names(bestmod), collapse = ", ")
  
  saveRDS(models, file = paste0(path, "sim_model_",j,"_", Sys.time(), ".rds"))
  
}

for (i in 1:100){
  test = as.vector(get(paste0("glint_iteration-",i)))
  if ("bmi, map, hdl, ltg, bmi:map, bmi:glu" %in% test){
    print(i)
  }
}
