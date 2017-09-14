dantzig<-function (X, Y, lambda = NULL, nlambda = NULL, lambda.min.value = NULL, 
  lambda.min.ratio = NULL, rho = 1, method = "lq", q = 2, 
  res.sd = FALSE, prec = 1e-05, max.ite = 1e+05, verbose = TRUE) 
{
  cat("Lu's dantzig\n")
  if (method != "dantzig" && method != "lq" && method != "lasso") {
    cat("\"method\" must be dantzig, lasso or lq.\n")
    return(NULL)
  }
  if (method == "lq") {
    if (q < 1 || q > 2) {
      cat("q must be in [1, 2] when method = \"lq\".\n")
      return(NULL)
    }
  }
  if (verbose) {
    cat("Sparse Linear Regression with L1 Regularization.\n")
  }
  if (sum(is.na(X)) > 0 || sum(is.na(Y)) > 0) {
    cat("The input has missing values.\n")
    Xrow.na = rowSums(is.na(X))
    nXrow.na = sum(Xrow.na)
    Xidx.na = which(Xrow.na > 0)
    Y.na = as.numeric(is.na(Y))
    nY.na = sum(Y.na)
    Yidx.na = which(Y.na > 0)
    idx.na = unique(c(Xidx.na, Yidx.na))
    X = X[-idx.na, ]
    Y = as.matrix(Y[-idx.na, ], ncol = 1)
    n = nrow(X)
    if (n == 0) {
      cat("Too many missing values.\n")
      return(NULL)
    }
  }
  if (is.null(X) || is.null(Y)) {
    cat("No data input.\n")
    return(NULL)
  }
  n = nrow(X)
  d = ncol(X)
  if (n == 0 || d == 0) {
    cat("No data input.\n")
    return(NULL)
  }
  maxdf = max(n, d)
  #xm = matrix(rep(colMeans(X), n), nrow = n, ncol = d, byrow = T)
  xm = matrix(0,n,d)
  
  x1 = X - xm
  #sdxinv = 1/sqrt(colSums(x1^2)/(n - 1))
  sdxinv = rep(1,d)
  xx = x1 * matrix(rep(sdxinv, n), nrow = n, ncol = d, byrow = T)
  #ym = mean(Y)
  ym = 0
  y1 = Y - ym
  if (res.sd == TRUE) {
    sdy = sqrt(sum(y1^2)/(n - 1))
    yy = y1/sdy
  }
  else {
    sdy = 1
    yy = y1
  }
  intercept = FALSE
  if (intercept) {
    xx = cbind(rep(1, nrow(xx)), xx)
    X = cbind(rep(1, nrow(X)), X)
    d = d + 1
  }
  if (!is.null(lambda)) 
    nlambda = length(lambda)
  if (is.null(lambda)) {
    if (is.null(nlambda)) 
      nlambda = 5
    if (method == "dantzig") {
      if (intercept) 
        lambda.max = max(abs(crossprod(xx[, 2:d], yy/n)))
      else lambda.max = max(abs(crossprod(xx, yy/n)))
    }
    if (method == "lq") {
      if (q == 2) {
        if (intercept) 
          lambda.max = max(abs(crossprod(xx[, 2:d], 
            yy/sqrt(sum(yy^2))/sqrt(n))))
        else lambda.max = max(abs(crossprod(xx, yy/sqrt(sum(yy^2))/sqrt(n))))
      }
      else {
        if (q == 1) {
          if (intercept) 
            lambda.max = max(abs(crossprod(xx[, 2:d], 
              sign(yy)/n)))
          else lambda.max = max(abs(crossprod(xx, sign(yy)/n)))
        }
        else {
          if (intercept) {
            lambda.max = max(abs(crossprod(xx[, 2:d], 
              sign(yy) * (abs(yy)^(q - 1))/(sum(abs(yy)^q)^((q - 
                1)/q))/n^(1/q))))
          }
          else {
            lambda.max = max(abs(crossprod(xx, sign(yy) * 
              (abs(yy)^(q - 1))/(sum(abs(yy)^q)^((q - 
              1)/q))/n^(1/q))))
          }
        }
      }
    }
    if (method == "lasso") {
      if (intercept) 
        lambda.max = max(abs(crossprod(xx[, 2:d], yy/n)))
      else lambda.max = max(abs(crossprod(xx, yy/n)))
    }
    if (method == "dantzig") {
      if (is.null(lambda.min.ratio)) {
        lambda.min.ratio = 0.5
      }
      if (is.null(lambda.min.value)) {
        lambda.min.value = lambda.min.ratio * lambda.max
      }
    }
    else {
      if (is.null(lambda.min.value)) {
        lambda.min.value = sqrt(log(d)/n)
      }
      else {
        if (is.null(lambda.min.ratio)) {
          lambda.min.ratio = lambda.min.value/lambda.max
        }
      }
    }
    if (lambda.max < lambda.min.value) {
      lambda.max = 1
      lambda.min.value = 0.4
    }
    lambda = exp(seq(log(lambda.max), log(lambda.min.value), 
      length = nlambda))
    rm(lambda.max, lambda.min.value, lambda.min.ratio)
    gc()
  }
  if (is.null(rho)) 
    rho = 1
  begt = Sys.time()
  if (method == "dantzig") {
    if (d >= n) 
      out = dantzig_inner(yy, xx, lambda, nlambda, 
        n, d, maxdf, rho, max.ite, prec, intercept, 
        verbose)
    else out = slim.dantzig.ladm.scr2(yy, xx, lambda, nlambda, 
      n, d, maxdf, rho, max.ite, prec, intercept, verbose)
    q = "infty"
  }
  if (method == "lq") {
    if (q == 1) 
      out = slim.lad.ladm.scr.btr(yy, xx, lambda, nlambda, 
        n, d, maxdf, rho, max.ite, prec, intercept, 
        verbose)
    if (q == 2) 
      out = slim.sqrt.ladm.scr(yy, xx, lambda, nlambda, 
        n, d, maxdf, rho, max.ite, prec, intercept, 
        verbose)
    if (q > 1 && q < 2) 
      out = slim.lq.ladm.scr.btr(yy, xx, q, lambda, nlambda, 
        n, d, maxdf, rho, max.ite, prec, intercept, 
        verbose)
  }
  if (method == "lasso") 
    out = slim.lasso.ladm.scr(yy, xx, lambda, nlambda, n, 
      d, maxdf, max.ite, prec, intercept, verbose)
  runt = Sys.time() - begt
  df = rep(0, nlambda)
  if (intercept) {
    for (i in 1:nlambda) df[i] = sum(out$beta[[i]][2:d] != 
      0)
  }
  else {
    for (i in 1:nlambda) df[i] = sum(out$beta[[i]] != 0)
  }
  est = list()
  intcpt0 = matrix(0, nrow = 1, ncol = nlambda)
  intcpt = matrix(0, nrow = 1, ncol = nlambda)
  if (intercept) {
    beta1 = matrix(0, nrow = d - 1, ncol = nlambda)
    for (k in 1:nlambda) {
      tmp.beta = out$beta[[k]][2:d]
      beta1[, k] = sdxinv * tmp.beta * sdy
      intcpt[k] = ym - as.numeric(xm[1, ] %*% beta1[, 
        k]) + out$beta[[k]][1] * sdy
      intcpt0[k] = intcpt[k]
    }
  }
  else {
    beta1 = matrix(0, nrow = d, ncol = nlambda)
    for (k in 1:nlambda) {
      tmp.beta = out$beta[[k]]
      intcpt0[k] = 0
      beta1[, k] = sdxinv * tmp.beta * sdy
      intcpt[k] = ym - as.numeric(xm[1, ] %*% beta1[, 
        k])
    }
  }
  est$beta0 = out$beta
  est$beta = beta1
  est$intercept = intcpt
  est$intercept0 = intcpt0
  est$Y = Y
  est$X = X
  est$lambda = lambda
  est$nlambda = nlambda
  est$df = df
  est$method = method
  est$q = q
  est$ite = out$ite
  est$verbose = verbose
  est$runtime = runt
  class(est) = "slim"
  if (verbose) 
    print(est)
  return(est)
}
