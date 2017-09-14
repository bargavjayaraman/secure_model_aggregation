my_sugm = function (data, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, 
          rho = NULL, method = "tiger", sym = "or", shrink = NULL, 
          prec = 1e-04, max.ite = 10000, standardize = FALSE, perturb = TRUE, 
          verbose = TRUE) 
{
  if (verbose) {
    cat("High-deimensional Sparse Undirected Graphical Models.\n")
  }
  if (method != "clime" && method != "tiger") {
    cat("\"method\" must be either \"clime\" or \"tiger\" \n")
    return(NULL)
  }
  n = nrow(data)
  d = ncol(data)
  if (method == "tiger" && d < 3) {
    cat("d>=3 is required for \"tiger\" \n")
    cat("More on help(sugm) \n")
    return(NULL)
  }
  if (method == "clime" && d < 2) {
    cat("d>=2 is required for \"clime\" \n")
    cat("More on help(sugm) \n")
    return(NULL)
  }
  maxdf = max(n, d)
  est = list()
  est$cov.input = isSymmetric(data)
  if (est$cov.input) {
    if (verbose) {
      cat("The input is identified as the covriance matrix.\n")
    }
    if (sum(is.na(data)) > 0) {
      cat("The input has missing values for covariance/correlation input.\n")
      return(NULL)
    }
    if (method == "tiger") {
      cat("The input for \"tiger\" cannot be covriance matrix.\n")
      return(NULL)
    }
    if (standardize) {
      S0 = data
      diag.cov = diag(S0)
      diag.cov.invsq = diag(1/sqrt(diag.cov))
      S = diag.cov.invsq %*% S0 %*% diag.cov.invsq
    }
    else {
      S0 = data
      S = S0
    }
  }
  if (!est$cov.input) {
    X0 = data
    if (method == "tiger" && sum(is.na(X0)) > 0) {
      cat("The input for \"tiger\" has missing values.\n")
      return(NULL)
    }
    X1 = X0 - matrix(rep(colMeans(X0), n), nrow = n, byrow = TRUE)
    S0 = crossprod(X1)/(n - 1)
    diag.cov = diag(S0)
    diag.cov.invsq = diag(1/sqrt(diag.cov))
    if (method == "tiger") {
      data = X1 %*% diag.cov.invsq
      S = S0
    }
    else {
      if (standardize) {
        S = cor(X0, use = "pairwise.complete.obs")
        diag.cov = diag(S)
        diag.cov.invsq = diag(1/sqrt(diag.cov))
      }
      else {
        S = cov(X0, use = "pairwise.complete.obs")
      }
    }
  }
  if (!is.null(lambda)) 
    nlambda = length(lambda)
  if (is.null(lambda)) {
    if (method == "tiger") {
      if (is.null(nlambda)) {
        nlambda = 5
      }
      if (is.null(lambda.min.ratio)) 
        lambda.min.ratio = 0.4
      lambda.max = pi * sqrt(log(d)/n)
      lambda = seq(lambda.max, lambda.min.ratio * lambda.max, 
                   length = nlambda)
    }
    else {
      if (is.null(nlambda)) 
        nlambda = 5
      if (is.null(lambda.min.ratio)) 
        lambda.min.ratio = 0.4
      lambda.max.tmp1 = min(max(S - diag(diag(S))), -min(S - 
                                                           diag(diag(S))))
      lambda.max.tmp2 = max(max(S - diag(diag(S))), -min(S - 
                                                           diag(diag(S))))
      if (lambda.max.tmp1 == 0) 
        lambda.max = lambda.max.tmp2
      else lambda.max = lambda.max.tmp1
      lambda.min = lambda.min.ratio * lambda.max
      lambda = exp(seq(log(lambda.max), log(lambda.min), 
                       length = nlambda))
      rm(lambda.max, lambda.min, lambda.min.ratio)
      gc()
    }
  }
  est$lambda = lambda
  est$nlambda = nlambda
  if (is.null(rho)) 
    rho = 1
  begt = Sys.time()
  if (method == "clime") {
    if (is.logical(perturb)) {
      if (perturb) {
        perturb = 1/sqrt(n)
        #perturb = 1e-5
      }
      else {
        perturb = 0
      }
    }
    S = S + diag(d) * perturb
    if (method == "clime") {
      if (is.null(shrink)) 
        shrink = 0
      if (is.null(max.ite)) 
        max.ite = 10000
      re.sugm = sugm.clime.ladm.scr(S, lambda, nlambda, 
                                    n, d, maxdf, rho, shrink, prec, max.ite, verbose)

    }
  }

  est$ite = re.sugm$ite
  runt = Sys.time() - begt
  for (j in 1:d) {
    if (re.sugm$col.cnz[j + 1] > re.sugm$col.cnz[j]) {
      idx.tmp = (re.sugm$col.cnz[j] + 1):re.sugm$col.cnz[j + 
                                                           1]
      ord = order(re.sugm$row.idx[idx.tmp])
      re.sugm$row.idx[idx.tmp] = re.sugm$row.idx[ord + 
                                                   re.sugm$col.cnz[j]]
      re.sugm$x[idx.tmp] = re.sugm$x[ord + re.sugm$col.cnz[j]]
    }
  }

  est$runtime = runt

  est$icov1 = re.sugm$icov1

  class(est) = "sugm"
  gc()
  return(est)
}
