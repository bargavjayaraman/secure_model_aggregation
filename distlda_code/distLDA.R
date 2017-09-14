# distributed LDA
# input data 2 * m list
distLDA <- function(data, lambda1, lambda2, t, isInDetail = FALSE, doTruncateForNaive = FALSE){
  #n = N / m
  #n2 = n - n1
  m = length(data) / 2
  d = nrow(data[[1]])
  betatilde = matrix(0, d, m)
  betahat = matrix(0, d, m)
  mu1 = matrix(0,d,m)
  mu2 = matrix(0,d,m)
  for (i in 1:m){
    tmp = SLDA(data[[i]], data[[i + m]], lambda1)
    mu1[, i] = tmp[[2]]
    mu2[, i] = tmp[[3]]
    betahat[, i] = tmp[[1]]
    tmp[[4]] = tmp[[4]] + diag(d) * 1e-3 # manual pertubation
    out = my_sugm(tmp[[4]], method = "clime", lambda = lambda2, perturb = FALSE)
    Theta = out$icov1[[1]]
    #out = fastclime(tmp[[4]], lambda.min = lambda2, nlambda = 100)
    #Theta = out$icovlist[[out$maxnlambda]]
    #if( max( out$lambdamtx[out$maxnlambda, ] ) > lambda2)
    #  cat("too small lambda2\n")
    #Ihat = crossprod(Theta,tmp[[4]])
    betatilde[, i] = betahat[, i] - crossprod(Theta, (tmp[[4]] %*% betahat[, i] - tmp[[2]] + tmp[[3]]))
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
