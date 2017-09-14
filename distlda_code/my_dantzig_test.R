## load library
library(flare)
source("my_dantzig.R")
source("my_dantzig_inner.R")
source("dantzig.R")
source("dantzig_inner.R")

## generate data
n = 100
d = 300
s = 30
X = matrix(rnorm(n*d), n, d)
beta = runif(d,0,1);
mask = (s+1):d
beta[mask]=0
eps = rnorm(n)
Y = X%*%beta + 1*eps
XTX = crossprod(X,X)/n
XTY = crossprod(X,Y)/n
nlamb = 5
ratio = 0.3
## Regression with "dantzig", general "lq" and "lasso" respectively
out1 = dantzig(X=X,Y=Y,lambda = 2*sqrt(log(d)/n),method="dantzig")
#out2 = slim(X=X,Y=Y,nlambda=nlamb,lambda.min.ratio=ratio,method="dantzig")
out2 = my_dantzig(XTX,XTY,lambda = 2*sqrt(log(d)/n),method="dantzig")
diff=out1$beta-out2$beta
print(norm(diff,'F'))
s1=which(out1$beta!=0)
s2=which(out2$beta!=0)
print(2.0*(length(intersect(s1,s2))+0.0)/(as.double(length(s1))+as.double(length(s2))))