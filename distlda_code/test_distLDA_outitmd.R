# test distLDA
# output intermediate result
library(flare)
#load("~/lab2T/dataset/syndata_10000_50.Rdata")
load("~/syndata_10000_100.Rdata")
n = N / m
beta = Theta %*% (mu1 - mu2)
expr_num = length(experiments)

source("distLDA.R")
source("centrLDA2.R")
source("SLDA.R")
source("my_sugm.R")
source("my_dantzig.R")
source("my_dantzig_inner.R")

clusterExport(cl, c("SLDA", "my_dantzig", "my_dantzig_inner", "my_sugm"))
worker.init = function(arg) {
  for(a in arg) library(a, character.only=TRUE)
  NULL
}
C1 = 3
C2 = 1
clusterCall(cl,worker.init,c("flare"))
tmps = parLapply(cl, experiments, distLDA, C1*sqrt(log(d)/n), C2*sqrt(log(d)/n), t=0,isInDetail = TRUE)
cat("distLDA\n")

C4 = 0.06
C5 = 0.0006

t = C4 * 9 * sqrt(log(d)/N) * 6 + C5 * 11 * 9^2 * log(d) * 6 / n # M=9,max(s,s')=11,|beta|_1=6

Sstar = which(beta!=0)
for (j in 1:expr_num) {
  barbeta = rowMeans(tmps[[j]][[1]])
  barbetahat = rowMeans(tmps[[j]][[2]])
  
  
  barbeta[abs(barbeta) <= t ] = 0
  barbetahat[abs(barbetahat) <= t ] = 0
}

