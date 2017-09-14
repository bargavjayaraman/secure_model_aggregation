my_dantzig<-function (X, Y, n,lambda = NULL, nlambda = NULL, lambda.min.value = NULL, 
  lambda.min.ratio = NULL, rho = 1, method = "lq", q = 2, 
  res.sd = FALSE, prec = 1e-05, max.ite = 1e+05, verbose = TRUE) #X is square, Y is d-dimensional
{
  cat("My dantzig\n")
  if (method != "dantzig" && method != "lq" && method != "lasso") {
    cat("\"method\" must be dantzig, lasso or lq.\n")
    return(NULL)
  }
  if (verbose) {
    cat("Sparse Linear Regression with L1 Regularization.\n")
  }
  if (sum(is.na(X)) > 0 || sum(is.na(Y)) > 0) {##with missing data
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
  
  d = ncol(X)
  if (n > d)
    n = d
  if (n == 0 || d == 0) {
    cat("No data input.\n")
    return(NULL)
  }
  maxdf = max(n, d)
  xx = X
  yy = Y
  intercept = FALSE
  if (!is.null(lambda)) 
    nlambda = length(lambda)
  if (is.null(lambda)) {
    if (is.null(nlambda)) 
      nlambda = 5
    if (method == "dantzig") 
      lambda.max = max(abs(Y))
    if (method == "dantzig") {
      if (is.null(lambda.min.ratio)) {
        lambda.min.ratio = 0.5
      }
      if (is.null(lambda.min.value)) {
        lambda.min.value = lambda.min.ratio * lambda.max
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
      out = my_dantzig_inner(yy, xx, lambda, nlambda, 
        n, d, maxdf, rho, max.ite, prec, intercept, 
        verbose)
    else out = my_dantzig_inner2(yy, xx, lambda, nlambda, 
      n, d, maxdf, rho, max.ite, prec, intercept, 
      verbose)
    q = "infty"
  }
  runt = Sys.time() - begt
  df = rep(0, nlambda)
  for (i in 1:nlambda) df[i] = sum(out$beta[[i]] != 0)

  est = list()
  intcpt0 = matrix(0, nrow = 1, ncol = nlambda)
  intcpt = matrix(0, nrow = 1, ncol = nlambda)
  beta1 = matrix(0, nrow = d, ncol = nlambda)
  for (k in 1:nlambda) {
      tmp.beta = out$beta[[k]]
      intcpt0[k] = 0
      beta1[, k] = tmp.beta
  }

  est$beta0 = out$beta
  est$beta = beta1
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
