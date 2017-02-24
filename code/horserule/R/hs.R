#' Horseshoe regression Gibbs-sampler
#'
#' Generates posterior samples using the horseshoe prior.
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
#' @importFrom mvnfast rmvn
hs = function(X,  y,  niter, hsplus=F, prior = NULL, thin=1, restricted= 0.001){
  #fast sampler
  #initialize non informative priors
  p = dim(X)[2]
  n = dim(X)[1]
  if (is.null(prior)){
  prior = rep(1, times=p)
  }
  X = as.matrix(X)
  sigma2 = 1
  lambda2 = runif(p)
  tau2 = 1
  nu = rep(1, times = p)
  eta2 = rep(1, times=p)
  psi = rep(1, times=p)
  xi = 1
  XtX = t(X)%*%X
  ###Fast sampler for big p
  #References:
  #Fast sampling with Gaussian scale-mixture priors in high-dimensional regression
  #A. Bhattacharya, A. Chakraborty, B. K. Mallick
  if(p>=n) {
    fastb = function(M, sigma2, Astar, y){
      p = dim(M)[2]
      phi = M/sqrt(sigma2)
      D = sigma2*Astar
      alpha = y/sqrt(sigma2)
      u = as.vector(rmvn(1, rep(0, times=p), D))
      delta = as.vector(rmvn(1, rep(0, times=n), diag(n)))
      v = phi%*%u + delta
      w = solve((phi%*%D%*%t(phi)+diag(n)), (alpha-v))
      u + D%*%t(phi)%*%w
    }
  } else {
  #####Rue Algorithm
  fastb = function(M, sigma2, Astar, y) {
      p = dim(M)[2]
      phi = M/sqrt(sigma2)
      alpha=y/sqrt(sigma2)
      D = sigma2*Astar
      diag(D) = 1/diag(D)
      L = t(chol(t(phi)%*%phi + D))
      v = solve(L, t(phi)%*%alpha)
      m = solve(t(L), v)
      z = as.vector(rmvn(1, rep(0, times=p), diag(p)))
      w = solve(t(L), z)
      m + w
    }
  }

  beta = matrix(0, nrow = niter/thin, ncol = p)
  lambda = matrix(0, nrow = niter/thin, ncol = p)
  sigma = c()
  tau = c()
  iter = 1
  samp = 1
  update = 1:p
  Lambda_star = diag(tau2*lambda2)
  b = fastb(X, sigma2, Lambda_star, y)
  scalelast = rep(1, times=p)
  for(iter in 1:niter){

    #sample beta
    if (restricted>0){
      scale = tau2*lambda2
      update = which((scale>restricted)|(scalelast>restricted))
      scalelast = scale
      Lambda_star = diag(scale[update])
      Xup = X[,update]
      b[update] = fastb(Xup, sigma2, Lambda_star, y)
    } else {

      Lambda_star = diag(tau2*lambda2)
      b = fastb(X, sigma2, Lambda_star, y)

    }

    #sample sigma2
    e = y - X%*%b
    shape = (n+p)/2
    scale = (t(e)%*%e)/2 + sum((b^2)/(tau2*lambda2))/2
    sigma2 = 1/rgamma(1,shape = shape, scale = 1/scale)

    #sample tau2
    shape = (p + 1)/2
    scale = 1/xi + sum((b^2)/lambda2)/(2*sigma2)
    tau2 = 1/rgamma(1, shape = shape, scale = 1/scale)

    #sample xi
    scale = 1 + 1/tau2
    xi = 1/rexp(1, scale)

    if (hsplus==T){
      #sample lambda2
      lambda2 = c()
      scale = 1/nu + (b^2)/(2*tau2*sigma2)
      for(l in 1:p){
        lambda2[l] = 1/rexp(1, scale[l])
      }

      #sample nu
      nu = c()
      scale = 1/eta2 + 1/lambda2
      for(l in 1:p){
        nu[l] = 1/rexp(1, scale[l])
      }


      #sample eta Horseshoe+
      eta2 = c()
      scale = 1/psi + 1/nu
      for(l in 1:p) {
        eta2[l] = 1/rexp(1, rate=scale[l])
      }
      #sample psi Horseshoe+
      psi = c()
      scale = 1/(prior^2) + 1/eta2
      for(l in 1:p){
        psi[l] = 1/rexp(1, rate = scale[l])
      }
    } else {
      lambda2 = c()
      scale = 1/nu + (b^2)/(2*tau2*sigma2)
      for(l in 1:p){
        lambda2[l] = 1/rexp(1, scale[l])
      }

      #sample nu
      nu = c()
      scale = 1/(prior^2) + 1/lambda2
      for(l in 1:p){
        nu[l] = 1/rexp(1, scale[l])
      }
    }




    ##store parameters
    ##thinning
    if (iter %% thin == 0){
      beta[samp,] = b
      sigma[samp] = sigma2
      tau[samp] = tau2
      lambda[samp,] = lambda2
      samp = samp+1
    }

    ##for diagnostics
    if( iter %% 1000 == 0){
      cat(paste("updated coefficients:", length(update)))
      cat(paste("iteration", iter, "complete.\n"))
    }
  }
  list(beta, sigma, tau, lambda)
}
