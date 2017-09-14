SLDA<-function(data1,data2,lambda) {
  d = nrow(data1)
  mu1hat = rowMeans(data1)
  n1 = ncol(data1)
  mu2hat = rowMeans(data2)
  n2 = ncol(data2)
  cdata1 = data1 - matrix(rep(mu1hat,n1),d,n1)
  cdata2 = data2 - matrix(rep(mu2hat,n2),d,n2)
  Sigmahat = (tcrossprod(cdata1,cdata1) + tcrossprod(cdata2,cdata2)) / (n1 + n2)
  out = my_dantzig(Sigmahat,mu1hat-mu2hat,n1+n2,lambda,method="dantzig")
  return (list(out$beta,mu1hat,mu2hat,Sigmahat))
}