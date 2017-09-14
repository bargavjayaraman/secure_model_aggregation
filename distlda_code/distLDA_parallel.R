# distributed LDA, parallelized implementation
# input data 2 * m list
library(doParallel)
distLDA_parallel <- function(data, lambda1, lambda2, t, isInDetail = FALSE, doTruncateForNaive = FALSE){
  m = length(data) / 2
  d = nrow(data[[1]])
  betatilde = matrix(0, d, m)
  betahat = matrix(0, d, m)
  mu1 = matrix(0,d,m)
  mu2 = matrix(0,d,m)
  out <- foreach(i=1:m, .combine=cbind) %dopar% {
    tmp = SLDA(data[[i]], data[[i + m]], lambda1)
    tmp[[4]] = tmp[[4]] + diag(d) * 1e-3 # manual pertubation
    out = my_sugm(tmp[[4]], method = "clime", lambda = lambda2, perturb = FALSE)
    Theta = out$icov1[[1]]
    Ihat = crossprod(Theta,tmp[[4]])
    bt = tmp[[1]] - crossprod(Theta, (tmp[[4]] %*% tmp[[1]] - tmp[[2]] + tmp[[3]]))/diag(Ihat)
    list(bt = bt, bh = tmp[[1]], m1 = tmp[[2]], m2 = tmp[[3]])
  }
  for(i in 1:m) {
    betatilde[, i] = out[[1, i]]
    betahat[, i] = out[[2, i]]
    mu1[, i] = out[[3, i]]
    mu2[, i] = out[[4, i]]
  }
  
  if(isInDetail)
    return(list(betatilde,betahat))
  else {
    barbeta = rowMeans(betatilde)
    barbetahat = rowMeans(betahat)
    barbeta[abs(barbeta) <= t ] = 0
    if(doTruncateForNaive) {
      barbetahat[abs(barbetahat) <= t] = 0
    }
    return(list(barbeta,barbetahat,rowMeans(mu1),rowMeans(mu2)))
  }
}
