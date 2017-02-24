#' predict.hs
#'
#' Predict unseen data with the horseshoe model.

#'@param object list containing a model of class "hs_rulefit".
#'@param newdata Dataframe containing the unseen data to predict.
#'@param burnin Number of samples that is disregarded as burnin. Increase number in case of slow convergence.
#'@param postmean If true returns the Predictive-Posterior mean value. If False returns the full predictive posterior distribution.
#'@param ... additional arguments
#' @return Returns the predictive posterior distribution matrix. Column i contains the predictive posterior of observation i.
#' @export

predict.HorseRulemodel = function(object, newdata, burnin=100, postmean=TRUE, ...) {
  linterms = object$modelstuff$linterms
  sdr = object$modelstuff$sdr
  mur = object$modelstuff$mur
  sdl = object$modelstuff$sdl
  mul = object$modelstuff$mul
  y = object$y
  ytransform = object$modelstuff$ytransform
  muy = object$modelstuff$muy
  sdy = object$modelstuff$sdy
  intercept = object$modelstuff$intercept

  beta = sdy*object$postdraws$beta[-c(1:burnin),]

  rules = object$rules

  Xrt = createXtest(newdata, rules)
  for(l in 1:dim(Xrt)[2]){
    Xrt[,l] = (Xrt[,l]-mur[l])/sdr[l]
  }

  ##preparing test data set. Standardize linear terms.
  if(length(linterms > 0)){
    for(l in 1:length(linterms)){
      newdata[,linterms[l]] = (newdata[,linterms[l]]-mul[l])/sdl[l]
    }
  }
  if(intercept==T){
    X_test = cbind(rep(1, times=dim(newdata)[1]), newdata[,linterms], Xrt)
  } else {
    X_test = cbind(newdata[,linterms], Xrt)
  }
  if(ytransform == "log") {
    predDist = apply(beta, 1, function(x)exp(as.matrix(X_test)%*%x + muy))
  } else if(is.numeric(y)) {
    predDist = apply(beta, 1, function(x)as.matrix(X_test)%*%x + muy)
  } else {
    predDist = apply(beta, 1, function(x)1/(1+exp(-(as.matrix(X_test)%*%x))))
  }
  if(postmean == T){
  apply(t(predDist),2,mean)
  } else {
  t(predDist)
  }

}
