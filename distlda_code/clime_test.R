## load package required
library(flare)
library(fastclime)
source("my_sugm.R")
## generating data
n = 25
d = 50
D = sugm.generator(n=n,d=d,graph="random",g=1)
hSigma = crossprod(D$data,D$data)/n
## sparse precision matrix estimation with method "clime"
out1 = sugm(hSigma, method = "clime",lambda = 0.3, perturb = FALSE)
out2 = my_sugm(hSigma, method = "clime",lambda = 0.3, perturb = FALSE)
out3 = fastclime(hSigma, lambda.min = 0.3, nlambda = 100)