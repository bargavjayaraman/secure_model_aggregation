## load library
library(flare)
source("dantzig.R")
source("dantzig_inner.R")
## generate data
n = 100
d = 200
s = 5
X = 2 * matrix(rnorm(n*d), n, d)
beta = runif(d,0,1);
mask = 6:d
beta[mask]=0
eps = rnorm(n)
Y = X%*%beta + 0.01*eps
nlamb = 5
ratio = 0.3
## Regression with "dantzig", general "lq" and "lasso" respectively
out1 = dantzig(X=X,Y=Y,nlambda=nlamb,lambda.min.ratio=ratio,method="dantzig")
out2 = slim(X=X,Y=Y,nlambda=nlamb,lambda.min.ratio=ratio,method="dantzig")
diff=out1$beta-out2$beta
print(norm(diff))
s1=which(out1$beta!=0)
s2=which(out2$beta!=0)
print(all(s1==s2))