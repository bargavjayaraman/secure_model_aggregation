# generate distributed LDA data for time complexity evaluation
N = 200
m = 1
expr_total = 20
n = N / m
n_int = as.integer(n)
N_res = N - n_int * m
alpha = 0.5
n1 = round(n_int * alpha)
n2 = n_int - n1
n1p = round((n_int + 1) * alpha)
n2p = n_int + 1 - n1p
d = 200

rho = 0.8
inva = 1 /(1 - rho^2)
invb = - rho * inva
invc = (1 + rho^2) / (1 - rho^2)
Sigma = matrix(0,d,d)
for (k in 1:d) {
  for (l in 1:d) {
    Sigma[k,l] = rho ^ abs(k - l)
  }
}

eigdec_result = eigen(Sigma,T)
#Theta = eigdec_result$vectors %*% tcrossprod(diag(1/eigdec_result$values), eigdec_result$vectors)
Theta = matrix(0, d, d)
for (k in 1:d){
  if (k == 1){
    Theta[k,k] = inva
    Theta[k,k+1] = invb
  }
  else if(k == d){
    Theta[k,k] = inva
    Theta[k,k - 1] = invb
  }
  else{
    Theta[k,k] = invc
    Theta[k,k - 1] = invb
    Theta[k,k + 1] = invb
  }
}
sqrtSigma = eigdec_result$vectors %*% tcrossprod(diag(sqrt(eigdec_result$values)), eigdec_result$vectors)

mu1 = matrix(rep(0,d),d,1)
mu2 = matrix(c(rep(1,10),rep(0,d-10)),d,1)

beta = Theta %*% (mu1 - mu2)

#randomize
experiments = vector("list", expr_total)
for (expr_no in 1:expr_total){
  rData = matrix(rnorm(N*d), d, N)
  rData = sqrtSigma %*% rData
  
  data = vector("list",2 * m)
  counter = 0
  for (i in 1:m) {
    if(i <= N_res) {
      data[[i]] = rData[, counter + (1:n1p)] + matrix(rep(mu1,n1p),d,n1p)
      counter = counter + n1p
      data[[i + m]] = rData[, counter + (1:n2p)] + matrix(rep(mu2,n2p),d,n2p)
      counter = counter + n2p
    } else {
      data[[i]] = rData[, counter + (1:n1)] + matrix(rep(mu1,n1),d,n1)
      counter = counter + n1
      data[[i + m]] = rData[, counter + (1:n2)] + matrix(rep(mu2,n2),d,n2)
      counter = counter + n2
    }
  }
  experiments[[expr_no]] = data
}
save(experiments,N,m,d,alpha,Theta,Sigma,mu1,mu2,file="~/syndata_200_1.Rdata")
