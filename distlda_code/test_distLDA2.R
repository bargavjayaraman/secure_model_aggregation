# test distLDA
# fix n, vary m, N
library(flare)
#load("~/lab2T/dataset/syndata_10000_50.Rdata")
load("~/syndata_10000_100.Rdata")
n = N / m
beta = Theta %*% (mu1 - mu2)
expr_num = length(experiments)
m_s = c(10,20,30,40,50)

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
C3 = 3
clusterCall(cl,worker.init,c("flare"))
tmps = parLapply(cl, experiments, distLDA, C1*sqrt(log(d)/n), C2*sqrt(log(d)/n), t=0,isInDetail = TRUE)
cat("distLDA\n")
betas = parLapply(cl, experiments, centrLDA2, C3*sqrt(log(d)/m_s/n), m_s)
cat("centrLDA\n")

C4 = 0.06
C5 = 0.0006
F1 = array(0, c(length(m_s), 3, expr_num))
Fnorms = array(0, c(length(m_s), 3, expr_num))
Mnorms = array(0, c(length(m_s), 3, expr_num))
Sstar = which(beta!=0)
for (i in 1:length(m_s)) {
  for (j in 1:expr_num) {
    tmp = rowMeans(tmps[[j]][[1]][,1:m_s[i]])
    barbetahat = rowMeans(tmps[[j]][[2]][,1:m_s[i]])
    barbeta = tmp
    t = C4 * 9 * sqrt(log(d)/n/m_s[i]) * 6 + C5 * 11 * 9^2 * log(d) * 6 / n # M=9,max(s,s')=11,|beta|_1=6
    barbeta[abs(barbeta) <= t ] = 0
    barbetahat[abs(barbetahat) <= t ] = 0
    Shat1 = which(betas[[j]][[i]] != 0) #centered
    Shat2 = which(barbetahat != 0) #naive
    Shat3 = which(barbeta != 0) #distlda
    precision1 = length(intersect(Sstar, Shat1))/length(Shat1)
    precision2 = length(intersect(Sstar, Shat2))/length(Shat2)
    precision3 = length(intersect(Sstar, Shat3))/length(Shat3)
    recall1 = length(intersect(Sstar, Shat1))/length(Sstar)
    recall2 = length(intersect(Sstar, Shat2))/length(Sstar)
    recall3 = length(intersect(Sstar, Shat3))/length(Sstar)
    F1[i,1,j] = 2*precision1*recall1 / (precision1 + recall1)
    if(is.nan(F1[i,1,j])) {
      F1[i,1,j] = 0
    }
    F1[i,2,j] = 2*precision2*recall2 / (precision2 + recall2)
    if(is.nan(F1[i,2,j])) {
      F1[i,2,j] = 0
    }
    F1[i,3,j] = 2*precision3*recall3 / (precision3 + recall3)
    if(is.nan(F1[i,3,j])) {
      F1[i,3,j] = 0
    }
    
    Fnorms[i,1,j]=norm(betas[[j]][[i]] - beta, 'F')
    Fnorms[i,2,j]=norm(barbetahat - beta, 'F')
    Fnorms[i,3,j]=norm(barbeta - beta, 'F')
    
    Mnorms[i,1,j]=norm(betas[[j]][[i]] - beta, 'M')
    Mnorms[i,2,j]=norm(barbetahat - beta, 'M')
    Mnorms[i,3,j]=norm(barbeta - beta, 'M')
  }
}
F1mean = apply(F1,c(1,2),mean)
F1std = apply(F1,c(1,2),sd)
Fnormsmean = apply(Fnorms,c(1,2),mean)
Fnormsstd = apply(Fnorms,c(1,2),sd)
Mnormsmean = apply(Mnorms,c(1,2),mean)
Mnormsstd = apply(Mnorms,c(1,2),sd)

save(F1,Fnorms,Mnorms,C1,C2,C3,C4,C5,file = "results/measure_10000_50_020114_t.Rdata")
