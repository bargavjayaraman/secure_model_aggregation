# misclassification rate in synthetic dataset for disributed results
load("syndata.Rdata")
Nexpr = 20
misclsfctn = matrix(0,5,Nexpr) #misclassification rate
m_s = c(20,40,60,80,100)
load("../../aaa_09192120.Rdata")
load("../../centr20000_100.Rdata")
for(stage in 1:5){
for(expr_no in 1:Nexpr){
  #betahat = read.table(paste('../Results_20000/20 parties/',as.character(expr_no),'.txt', sep=''))
  #betahat = as.numeric(betahat)
  betahat = rowMeans(tmps[[expr_no]][[1]][, 1:m_s[stage]])
  #betahat = betas[[expr_no]][[stage]]
  
  # if using 40 parties, should be <=t[2], 60 parties for t[3], 80 for t[4] ...
  betahat[abs(betahat) <= t[stage]/2 ] = 0
  
  # if using 40 parties, should be muhat[[expr_no]][, 2], 60 for muhat[[expr_no]][, 3]...
  difference = apply(data[[1]], 2, '-', muhat[[expr_no]][, stage])
  
  
  
  predict = (crossprod(betahat, difference) > 0)
  e1 = sum(predict == F)
  
  # if using 40 parties, should be muhat[[expr_no]][, 2], 60 for muhat[[expr_no]][, 3]...
  difference = apply(data[[2]], 2, '-', muhat[[expr_no]][, stage])
  
  predict = (crossprod(betahat, difference) > 0)
  e2 = sum(predict == T)
  misclsfctn[stage,expr_no] = (e1 + e2) / (10000)
}
}
cat( rowMeans(misclsfctn) )
cat('\n')
cat( apply(misclsfctn,1,sd) )