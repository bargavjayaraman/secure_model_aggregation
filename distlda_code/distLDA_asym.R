# distributed LDA, assymetric
# different machines hold different number of data
# Therefore we need individual lambda1/lambda2 for each machine
# input data 2 * m list
distLDA_asym <- function(data, lambda1, lambda2, t, isInDetail = FALSE, doTruncateForNaive = FALSE){
  #n = N / m
  #n2 = n - n1
  m = length(data) / 2
  d = nrow(data[[1]])
  betatilde = matrix(0, d, m)
  betahat = matrix(0, d, m)
  mu1 = matrix(0,d,m)
  mu2 = matrix(0,d,m)
  n1_vec = rep(0,m)
  n2_vec = rep(0,m)
  for (i in 1:m){
    tmp = SLDA(data[[i]], data[[i + m]], lambda1[i])
    mu1[, i] = tmp[[2]]
    mu2[, i] = tmp[[3]]
    betahat[, i] = tmp[[1]]
    tmp[[4]] = tmp[[4]] + diag(d) * 1e-3 # manual pertubation
    out = my_sugm(tmp[[4]], method = "clime", lambda = lambda2[i], perturb = FALSE)
    Theta = out$icov1[[1]]
    n1_vec[i] = ncol(data[[i]])
    n2_vec[i] = ncol(data[[i+m]])
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
    barbeta = apply(betatilde,1,weighted.mean,n1_vec+n2_vec)
    barbetahat = apply(betahat,1,weighted.mean,n1_vec+n2_vec)
    barbeta[abs(barbeta) <= t ] = 0
    if(doTruncateForNaive) {
      barbetahat[abs(barbetahat) <= t] = 0
    }
    return(list(barbeta,barbetahat,apply(mu1,1,weighted.mean,n1_vec),
                apply(mu2,1,weighted.mean,n2_vec)))
  }
}
