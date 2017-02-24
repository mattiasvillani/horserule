#' convergence_plot
#'
#' Can be used to check model convergence.

#'@param model list containing a model of class "hs_rulefit".
#'@param Xtest Out of bag sample to check error.
#'@param ytest response of test data.
#'@param burnin Number of samples disregarded as burnin.
#'@details To check convergence on the Training data just use training X and y instead of Xtest and ytest.

#' @export
#' @import graphics

convergence_plot = function(model, Xtest, ytest, burnin = 0){
  predDist = predict(model, Xtest)
  tau = model[[2]][[3]]
  sigma = model[[2]][[2]]
  erroreach = c()
  errortotal = c()
  if(is.numeric(ytest)){
  for(i in 1:dim(predDist)[1]){
    erroreach[i] = sqrt(mean((predDist[i,]-ytest)^2))
    if(i>1){
      pred = apply(predDist[1:i,], 2, mean)
      errortotal[i] = sqrt(mean((pred-ytest)^2))
    } else {
      errortotal[i] = erroreach[i]
    }
  }
  } else {
    ytest = as.numeric(ytest)-1
  for(i in 1:dim(predDist)[1]){
      erroreach[i] = 1-mean(ifelse(predDist[i,]<0.5,0,1) == ytest)
      if(i>1){
        pred = apply(predDist[1:i,], 2, mean)
        errortotal[i] = 1-mean(ifelse(pred<0.5,0,1)==ytest)
    }else {
        errortotal[i] = 1-erroreach[i]
    }
  }
  }
  layout(matrix(1:4, ncol=2))
  if(burnin == 0){
  plot(erroreach, type = "l", xlab = "Iteration", ylab = "RMSE",ylim=c(min(min(erroreach),min(errortotal)),max(max(erroreach),max(errortotal))))
  lines(errortotal, col = "orange", lwd = 3)
  plot(tau, type = "l", xlab="Iteration", ylab = "tau")
  if(is.numeric(ytest)){
  plot(sigma, type = "l", xlab="Iteration", ylab = "sigma")
  }
  } else {
  plot(erroreach[-c(1:burnin)], type = "l", xlab = "Iteration", ylab = "RMSE",ylim=c(min(min(erroreach),min(errortotal)),max(max(erroreach),max(errortotal))))
  lines(errortotal[-c(1:burnin)], col = "orange", lwd = 3)
  plot(tau[-c(1:burnin)], type = "l", xlab="Iteration", ylab = "tau")
  if(is.numeric(ytest)){
  plot(sigma[-c(1:burnin)], type = "l", xlab="Iteration", ylab = "sigma")
  }
  }
}
