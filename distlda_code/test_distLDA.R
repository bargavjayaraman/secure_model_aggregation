# test distLDA
library(flare)
load("~/lab2T/dataset/syndata_10000_100.Rdata")

source("distLDA.R")
source("centrLDA.R")
source("SLDA.R")
source("my_sugm.R")
source("my_dantzig.R")
source("my_dantzig_inner.R")
n = N/m

C1 = 3
C2 = 1
C3 = 3

beta = Theta %*% (mu1 - mu2)
expr_num = length(experiments)
Fnorms = matrix(0,3,expr_num)
Mnorms = matrix(0,3,expr_num)
F1 = matrix(0,3,expr_num)
#t = C4 * 9 * sqrt(log(d)/N) * 6 + C5 * 11 * 9^2 * m * log(d) * 6 / N
t = 0.1

#cl = makeCluster(detectCores())
clusterExport(cl, c("SLDA", "my_dantzig", "my_dantzig_inner", "my_sugm"))
worker.init = function(arg) {
  for(a in arg) library(a, character.only=TRUE)
  NULL
}
clusterCall(cl,worker.init,c("flare"))

beta3s = parLapply(cl, experiments, centrLDA, C3*sqrt(log(d)/N))
cat("centrLDA\n")
tmps = parLapply(cl, experiments, distLDA, C1*sqrt(log(d)/n), C2*sqrt(log(d)/n), t, doTruncateForNaive = TRUE)
cat("distLDA\n")

##########
Sstar = which(beta!=0)

for (expr_no in 1:expr_num){
  Shat1 = which(beta3s[[expr_no]] != 0) #centralized
  Shat2 = which(tmps[[expr_no]][[2]]!=0) #naive
  Shat3 = which(tmps[[expr_no]][[1]]!=0) #debiased
  precision1 = length(intersect(Sstar, Shat1))/length(Shat1)
  precision2 = length(intersect(Sstar, Shat2))/length(Shat2)
  precision3 = length(intersect(Sstar, Shat3))/length(Shat3)
  recall1 = length(intersect(Sstar, Shat1))/length(Sstar)
  recall2 = length(intersect(Sstar, Shat2))/length(Sstar)
  recall3 = length(intersect(Sstar, Shat3))/length(Sstar)
  F1[1,expr_no] = 2*precision1*recall1 / (precision1 + recall1)
  F1[2,expr_no] = 2*precision2*recall2 / (precision2 + recall2)
  F1[3,expr_no] = 2*precision3*recall3 / (precision3 + recall3)
  if(is.nan(F1[1,expr_no])) {
    F1[1,expr_no] = 0
  }
  if(is.nan(F1[2,expr_no])) {
    F1[2,expr_no] = 0
  }
  if(is.nan(F1[3,expr_no])) {
    F1[3,expr_no] = 0
  }
  Fnorms[3,expr_no]=norm(tmps[[expr_no]][[1]] - beta, 'F')
  Fnorms[2,expr_no]=norm(tmps[[expr_no]][[2]] - beta, 'F')
  Fnorms[1,expr_no]=norm(beta3s[[expr_no]] - beta, 'F')

  Mnorms[3,expr_no]=norm(tmps[[expr_no]][[1]] - beta, 'M')
  Mnorms[2,expr_no]=norm(tmps[[expr_no]][[2]] - beta, 'M')
  Mnorms[1,expr_no]=norm(beta3s[[expr_no]] - beta, 'M')
}

F1mean = rowMeans(F1)
F1std = apply(F1,1,sd)

Fnormsmean = rowMeans(Fnorms)
Fnormsstd = apply(Fnorms,1,sd)

Mnormsmean = rowMeans(Mnorms)
Mnormsstd = apply(Mnorms,1,sd)

save(tmps, beta3s, F1, Fnorms, Mnorms,C1,C2,C3,t,file="results/criteria_10000_100_020223.Rdata")
