#' Horseshoe regression Gibbs-sampler for classification
#'
#' Generates posterior samples using the horseshoe prior using the polya gamma data augmentation approach for logistic regression.
#' @param X A matrix containing the predictor variables to be used.
#' @param y The vector of numeric responses.
#' @param niter Number of posterior samples.
#' @param hsplus If "hsplus=T" the horseshoe+ extension will be used.
#' @param prior Prior for the individual predictors. If all 1 a standard horseshoe model is fit.
#' @param thin If > 1 thinning is performed to reduce autocorrelation.
#' @param restricted Threshold for restricted Gibbs sampling. In each iteration only coefficients with scale > restricted are updated. Set restricted = 0 for unrestricted Gibbs sampling.
#' @return A list containing the posterior samples of the following parameters:
#' \item{beta}{Matrix containing the posterior samples for the regression coefficients.}
#'  \item{sigma}{Vector contraining the Posterior samples of the error variance.}
#'  \item{tau}{Vector contraining the Posterior samples of the overall shrinkage.}
#' \item{lambda}{Matrix containing the posterior samples for the individual shrinkage parameter.}
#' @export
#' @import BayesLogit
#' @importFrom mvnfast rmvn
#'
hs_class = function(X,  y,  niter, hsplus=F, prior, thin=1, restricted= 0.001){
  p = dim(X)[2]
  n = dim(X)[1]
  X = as.matrix(X[,-1])
  sigma2 = 1
  lambda2 = runif(p-1)
  tau2 = 1
  nu = rep(1, times = p-1)
  xi = 1
  b = rep(0, times = p)

  ###Fast sampler for big p
  #References:
  #Fast sampling with Gaussian scale-mixture priors in high-dimensional regression
  #A. Bhattacharya, A. Chakraborty, B. K. Mallick
  fastb = function(X, sigma2, Astar, ystar){
    phi = X/sqrt(sigma2)
    D = sigma2*Astar
    alpha = ystar/sqrt(sigma2)
    u = as.vector(rmvn(1, rep(0, times=p), D))
    delta = as.vector(rmvn(1, rep(0, times=n), diag(n)))
    v = phi%*%u + delta
    w = solve((phi%*%D%*%t(phi)+diag(n)), (alpha-v))
    u + D%*%t(phi)%*%w
  }
  robustinv = function(X) {
    tryCatch(V<-solve(X),
             error = function(e) {print(paste("Inverse fail", 5)); as.integer(9)})
  }
  #storage matrices
  beta = matrix(0, nrow = niter/thin, ncol = p)
  lambda = matrix(0, nrow = niter/thin, ncol = p-1)
  sigma = c()
  tau = c()
  iter = 1
  samp = 1
  y = as.numeric(y) - 1

  n0 = sum(y==0)
  n1 = sum(y==1)

  z = vector(mode="numeric", length=n)
  for(iter in 1:niter){

    Lambda_star_inv = diag(1/(tau2*lambda2))
    Xb = b[1] + X%*%b[-1]
    w = sapply(Xb,function(x) rpg(1, 1, x))
    kappa = y-0.5
    z = kappa/w

    sig = diag(w)
    V = robustinv(t(X)%*%sig%*%X + Lambda_star_inv)
    if (is.matrix(V)){
      b[-1] = as.vector(rmvn(1, V%*%t(X)%*%sig%*%(z-b[1]*rep(1, n)), V))
    }
    v = (z - Xb)*w
    s = sum(w)

    b[1] = rnorm(1,(1/s)*sum(v), 1/s)


    lambda2 = c()
    scale = 1/nu + (b[-1]^2)/(2*tau2*sigma2)
    for(l in 1:(p-1)){
      lambda2[l] = 1/rexp(1, scale[l])
    }

    #sample tau2
    shape = (p+1)/2
    scale = 1/xi + sum((b[-1]^2)/lambda2)/(2*sigma2)
    tau2 = 1/rgamma(1, shape = shape, scale = 1/scale)

    #sample nu
    nu = c()
    scale = 1/(prior[-1]^2) + 1/lambda2
    for(l in 1:(p-1)){
      nu[l] = 1/rexp(1, scale[l])
    }
    #sample xi
    scale = 1 + 1/tau2
    xi = 1/rexp(1, scale)

    ##store parameters
    ##thinning
    if (iter %% thin == 0){
      beta[samp,] = b
      tau[samp] = tau2
      lambda[samp,] = lambda2
      samp = samp+1
    }
    ##for diagnostics
    if( iter %% 100 == 0)cat(paste("iteration", iter, "complete.\n"))
  }
  list(beta, sigma, tau, lambda)
}
