# misclassification rate in synthetic dataset for centralized
load("syndata.Rdata")
Nexpr = 20
misclsfctn = rep(0,Nexpr) #misclassification rate
load("../../centr20000_100.Rdata")
for(expr_no in 1:Nexpr){
  betahat = betas[[expr_no]][[5]]
  betahat = as.numeric(betahat)
  
  # if using 40 parties, should be <=t[2], 60 parties for t[3], 80 for t[4] ...
  #betahat[abs(betahat) <= t[1]/2 ] = 0
  
  # if using 40 parties, should be muhat[[expr_no]][, 2], 60 for muhat[[expr_no]][, 3]...
  difference = apply(data[[1]], 2, '-', muhat[[expr_no]][, 5])
  
  
  
  predict = (crossprod(betahat, difference) > 0)
  e1 = sum(predict == F)
  
  # if using 40 parties, should be muhat[[expr_no]][, 2], 60 for muhat[[expr_no]][, 3]...
  difference = apply(data[[2]], 2, '-', muhat[[expr_no]][, 5])
  
  predict = (crossprod(betahat, difference) > 0)
  e2 = sum(predict == T)
  misclsfctn[expr_no] = (e1 + e2) / (10000)
}

cat( mean(misclsfctn) )
cat('\n')
cat( sd(misclsfctn) )