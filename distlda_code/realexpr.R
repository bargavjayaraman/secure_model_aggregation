# experiments on real dataset
library(R.matlab)
library(flare)
source("distLDA_parallel.R")
source("distLDA.R")
source("SLDA.R")
source("my_sugm.R")
source("my_dantzig.R")
source("my_dantzig_inner.R")

worker.init = function(arg) {
  for(a in arg) library(a, character.only=TRUE)
  NULL
}
clusterCall(cl,worker.init,c("flare"))
clusterExport(cl, c("SLDA", "my_dantzig", "my_dantzig_inner", "my_sugm"))

a = readMat("~/lab2T/dataset/TDT2.mat")
class1idx = which(a$gnd==1)
class2idx = which(a$gnd==2)
N = 1200
traindata1 = sample(class1idx, N)
traindata2 = sample(class2idx, N)
testdata1 = setdiff(class1idx, traindata1)
testdata2 = setdiff(class2idx, traindata2)

m = 10
n = N/m
d = 200
n_int = floor(n)
N_res = N - m * n_int
expr = vector("list", 2*m)
counter = 0
for(i in 1:m) {
  if(i<=N_res) {
    expr[[i]] = as.matrix(t(a$fea[traindata1[counter+(1:(n_int+1))],1:d]))
    expr[[i+m]] = as.matrix(t(a$fea[traindata2[counter+(1:(n_int+1))],1:d]))
    counter = counter + n_int + 1
  }
  else {
    expr[[i]] = as.matrix(t(a$fea[traindata1[counter+(1:n_int)],1:d]))
    expr[[i+m]] = as.matrix(t(a$fea[traindata2[counter+(1:n_int)],1:d]))
    counter = counter + n_int
  }
}
para=list()
para$lambda1 = 3*sqrt(log(d)/n/2)
para$lambda2 = sqrt(log(d)/n/2)
para$t = 0.5
out1=distLDA_parallel(expr, para$lambda1, para$lambda2,para$t)
muhat = (out1[[3]]+out1[[4]]) / 2

testresult1 = ( out1[[1]] %*% ( t(a$fea[testdata1,1:d]) - matrix(rep(muhat,length(testdata1)),d,length(testdata1)) ) ) > 0

testresult2 = ( out1[[1]] %*% ( t(a$fea[testdata2,1:d]) - matrix(rep(muhat,length(testdata2)),d,length(testdata2)) ) ) > 0

testresult_naive1 = ( out1[[2]] %*% ( t(a$fea[testdata1,1:d]) - matrix(rep(muhat,length(testdata1)),d,length(testdata1)) ) ) > 0

testresult_naive2 = ( out1[[2]] %*% ( t(a$fea[testdata2,1:d]) - matrix(rep(muhat,length(testdata2)),d,length(testdata2)) ) ) > 0
