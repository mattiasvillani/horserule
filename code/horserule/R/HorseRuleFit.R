#' Horseshoe RuleFit
#'
#' fits a horseshoe rulefit model
#' @param model.formula formula type argument specifying the linear model. All standard lm arguments can be passed over, such as interactions and cubic terms.
#' @param data If model.forumla is used data must be the data frame containing the variables.
#' @param X A matrix or dataframe containing the predictor variables to be used.
#' @param y A vector containing the response variables. If numeric regression is performed and classification otherwise.
#' @param Xtest optional matrix or dataframe containing predictor variables of test set.
#' @param ytest optional vector containing the response values of the test set.
#' @param niter number of iterations for the horseshoe sampling.
#' @param burnin number of initial samples to be disregarded as burnin.
#' @param thin thinning parameter.
#' @param restricted Threshold for restricted Gibbs sampling. In each iteration only coefficients with scale > restricted are updated. Set restricted = 0 for unrestricted Gibbs sampling.
#' @param shrink.prior Specifies the shrinkage prior to be used for regularization. Currently the options "HS" and "HS+" for the Horseshoe+ are supported.
#' @param beta Hyperparameter to control the extra shrinkage on the rule complexity meassured as the rule length.
#' @param alpha Hyperparameter to control the extra shrinkage on the rules that cover only few observations. Set alpha = beta = 0 for the standard horseshoe without rule structure prior.
#' @param linp Hyperparameter to control prior shrinkage of linear terms. Set linp > 1 if strong linear effects are assumed.
#' @param ntree Number of trees in the ensemble step from which the rules are extracted.
#' @param ensemble Which ensemble method should be used to generate the rules? Options are "RF","GBM" or "both".
#' @param mix If ensemble = "both" mix*ntree are generated via random forest and (1-mix)*ntree trees via gradient boosting.
#' @param L Parameter controling the complexity of the generated rules. Higher values lead to more complex rules.
#' @param S Parameter controlling the minimum number of observations in the tree growing process.
#' @param minsup Rules with support < minsup are removed. Can be used to prevent overfitting.
#' @param linterms specifies the columns in X which should be included as linear terms in the hs rulefit model. Specified columns need to be numeric. Categorical variables have to be transformed (e.g. to dummies) before included as linear effects.
#' @param intercept If TRUE an intercept is included. Note that the y by default is shifted to have 0 mean therefor not necessary for regression. For classification highly recommended.
#' @param ytransform Choose "log" for logarithmic transform of y.
#' @return An object of class HorseRuleFit, which is a list of the following components:
##'
##'  \item{bhat}{Posterior mean of the regression coefficients.}
##'  \item{posteriorsamples}{List contraining the Posterior samples of the regression coefficients, error variance sigma and shrinkage tau.}
##'  \item{rules}{Vector containing the decision rules.}
##'  \item{Xt}{Matrix of train data with rules as dummies.}
##'  \item{y}{Response in train data.}
##'  \item{prior}{Vector rule structure prior for the individual rules.}
##'  \item{modelstuff}{List contraining the parameters used and values used for the normalization (means and sds).}
##'  \item{pred}{If Test data was supplied, gives back the predicted values.}
##'  \item{err}{If y-test was also supplies additionally gives back a test error score (RMSE for regression, Missclassificationrate for Classficitaion).}

#' @examples
#'library(MASS)
#'library(horserule)
#'data(Boston)
#' # Split in train and test data
#'N = nrow(Boston)
#'train = sample(1:N, 400)
#'Xtrain = Boston[train,-14]
#'ytrain = Boston[train, 14]
#'Xtest = Boston[-train, -14]
#'ytest = Boston[-train, 14]
#'
#' # Run the HorseRuleFit with 200 trees. Increase Number for better modelfit.
#' \dontrun{
#'hrres = HorseRuleFit(X = Xtrain, y=ytrain,
#'                     thin=1, niter=1000, burnin=100,
#'                     L=5, S=6, ensemble = "both", mix=0.3, ntree=200,
#'                     intercept=FALSE, linterms=1:13, ytransform = "log",
#'                     alpha=1, beta=2, linp = 1, restricted = 0)
#'
#' # Calculate the error
#'pred = predict(hrres, Xtest, burnin=100, postmean=TRUE)
#'sqrt(mean((pred-ytest)^2))
#'
#' # Look at the most important rules/linear effects.
#' importance_hs(hrres)
#'
#' # Look at the input variable importance.
#'Variable_importance(hrres, var_names=colnames(Xtrain))
#'}
#' @export
#' @import stats
#' @import utils
HorseRuleFit = function(model.formula=NULL,data=NULL,X=NULL, y=NULL, Xtest=NULL, ytest=NULL,
                        niter=1000,burnin=100, thin=1, restricted=0.001, shrink.prior ="HS",
                        beta=2, alpha=1, linp=1,
                        ensemble= "RF", L=4, S=6, ntree=250, minsup=.025, mix=0.5,
                        linterms=NULL, intercept=F, ytransform = "linear") {
  sdy=1
  muy=0
  inputtype = "matrix-type"
  terms_test = NULL
  terms_train = NULL
  if(is.null(model.formula)==F){
    if(is.null(data)){
      stop("If formula type argument is used, the argument 'data' needs to be specified.")
    } else{
      inputtype = "formula-type"
      yname = as.character(attr(terms(model.formula,data=data), "variables")[[attr(terms(medv~.,data=data), "response")+1]])
      terms_test = delete.response(terms(model.formula, data=data))
      attr(terms_test,"intercept") = 0
      if(is.null(yname)==F){
        eval(parse( text = (paste("y <- data$",yname, sep=""))))
        terms_train = terms(model.formula, data=data)
        attr(terms_train,"intercept") = 0
        X = model.matrix(terms_train, data=data)
      } else {
        stop("Error in the given formula, response not found in dataframe.")
      }

    }
  }

  if((is.matrix(X)|is.data.frame(X))==F){
    stop("X must be a matrix or data frame")
  }
  if((!is.numeric(y))&(length(unique(y))>2)){
    stop("y is not numeric and has more than 2 categories. Currently only regression and binary classification are supported.")
  }
  if(is.null(Xtest)==F){
    if(dim(X)[2]!=dim(Xtest)[2]){
      stop("The dimensionality between X and Xtest differ.")
    }
  }
  if(is.null(ytest)==F){
    if(mode(y)!=mode(ytest)){
      stop("The mode of y and ytest differs.")
    }
  }
  if(niter<=burnin){
    stop("Number of iterations needs to be higher than number of burnin.")
  }
  if(L<2){
    stop("Parameter L needs to be >=2.")
  }
  if(S<1){
    stop("ParameterS needs to be >=1.")
  }
  if(ntree<2){
    stop("Too few trees are chosen for Hs-RuleFit. For the standard horseshoe regression, please use hs() instead.")
  }
  if((minsup<0)|(minsup>=1)){
    stop("invalid choice for minimum support, please chose a value between 0 and 1.")
  }
  if((mix<0)|(mix>=1)){
    stop("invalid choice for mix, please chose a value between 0 and 1.")
  }

  if(is.null(linterms)){
    linear = F
  } else {
    linear = T
    if((is.numeric(linterms)==F) & (is.integer(linterms)==F)){
      stop("Invalid linterms. Must be a vector of type either numeric or integer, enumerating the linear terms to include.")
    }
    for(l in 1:length(linterms)){
      if(is.numeric(X[,linterms[l]])==F){
        stop(sprintf("Variable %i is not numeric and can not be included as linear term. Please check the variable.",l))
      }
    }
  }

  if(is.logical(intercept)==F){
    stop("Invalid intercept choice. Must be TRUE or FALSE.")
  }
  if(!(ytransform%in%c("linear", "log"))){
    stop("Not supported ytransform. Only logarithmic transform is supported currently.")
  }

  if(ytransform == "log") {
    y = log(y)
    if(any(is.na(y))){
      stop("logarithmic transform produced NAs. Please consider adding a small amount to the y before rerunning.")
    }
  }
  muy = 0

  if(is.numeric(y)){
    if(mean(y) != 0) {
      muy = mean(y)
      sdy = sd(y)
      yz = (y-muy)/sdy
    }
  } else {
    yz = y
  }

  N = length(y)
  if (ensemble == "RF") {
    capture.output(rulesf <- genrulesRF(X, yz, nt=ntree, S=S, L=L), file='NUL')
  } else if (ensemble == "GBM") {
    capture.output(rulesf <- genrulesGBM(X, yz, nt=ntree,S=S, L=L), file='NUL')
  } else if (ensemble == "both"){
    capture.output(rules1 <- genrulesRF(X, yz, nt=round(ntree*mix), S=S, L=L), file='NUL')
    capture.output(rules2 <- genrulesGBM(X, yz, nt=round(ntree*(1-mix)), S=S, L=L), file='NUL')
    rulesf = c(rules1, rules2)
  } else {
    print("invalid Tree ensemble choice")
    break
  }

  dt = createX(X= X, rules = rulesf, t = minsup)
  Xr = dt[[1]]
  rulesFin = dt[[2]]

  prior = calc_prior(rules=rulesFin, Xr,alpha=alpha, beta=beta)
  if(intercept == T){
    prior = c(1, rep(linp, times=length(linterms)), prior)
  } else {
    prior = c(rep(linp, times=length(linterms)), prior)
  }

  #standardize the rule frame
  mur = apply(Xr, 2, mean)
  sdr = apply(Xr, 2, sd)
  for(l in 1:dim(Xr)[2]){
    Xr[,l] = (Xr[,l]-mur[l])/sdr[l]
  }

  sdl=0
  mul=0

  if(length(linterms)>1){
    mul = apply(X[,linterms], 2, mean)
    sdl = apply(X[,linterms], 2, sd)
    for(l in 1:length(linterms)){
      X[,linterms[l]] = (X[,linterms[l]]-mul[l])/sdl[l]
    }
  } else if(length(linterms)==1) {
    mul = mean(X[,linterms])
    sdl = sd(X[,linterms])
    X[,linterms] = (X[,linterms] - mul)/sdl
  }



  if(linear==FALSE){
    if(intercept==TRUE){
      Xt = cbind(rep(1, times= dim(Xr)[1]),Xr)
    } else {
      Xt = Xr
    }
  } else {
    if(intercept==TRUE){
      Xt = cbind(rep(1, times=dim(X)[1]),X[,linterms], Xr)
    } else {
      Xt = cbind(X[,linterms], Xr)
    }
  }

  if (shrink.prior == "HS") {
    if(is.numeric(y)){
      hsmodel = hs(X = Xt, y=yz, niter=niter, prior=prior,thin=thin, hsplus=F, restricted=restricted)
      beta = sdy*hsmodel[[1]][-c(1:(burnin/thin)),]
      bhat = apply(beta, 2, mean)
    } else {
      hsmodel = hs_class(X = Xt, y=yz, niter=niter, prior=prior,thin=thin, hsplus=F, restricted=restricted)
      beta = sdy*hsmodel[[1]][-c(1:(burnin/thin)),]
      bhat = apply(beta, 2, mean)
    }
  } else if (shrink.prior =="HS+") {
    if(is.numeric(y)){
      hsmodel = hs(X = Xt, y=yz, niter=niter, prior=prior,thin=thin, hsplus=T, restricted=restricted)
      beta = sdy*hsmodel[[1]][-c(1:(burnin/thin)),]
      bhat = apply(beta, 2, mean)
    } else {
      hsmodel = hs_class(X = Xt, y=yz, niter=niter, prior=prior,thin=thin, hsplus=T, restricted=restricted)
      beta = sdy*hsmodel[[1]][-c(1:(burnin/thin)),]
      bhat = apply(beta, 2, mean)
    }
  } else {
    print("invalid prior choice")
    break
  }


  if(is.null(Xtest)==T){
    out = list(postmean=bhat,postdraws=list(beta=hsmodel[[1]], sigma=hsmodel[[2]], tau=hsmodel[[3]]), rules=rulesFin, df=Xt, y=y, prior=prior, modelstuff=list(linterms=linterms, sdr=sdr, mur=mur, sdl=sdl, mul=mul, ytransform=ytransform, muy=muy,sdy =sdy,sdl=sdl, intercept=intercept, inputtype=inputtype, model.formula = model.formula, terms_test = terms_test), X=X)
    class(out) = "HorseRulemodel"
  } else {
    ##create rules. Standardize.
    Xrt = createXtest(Xtest, rulesFin)
    for(l in 1:dim(Xr)[2]){
      Xrt[,l] = (Xrt[,l]-mur[l])/sdr[l]
    }

    ##preparing test data set. Standardize linear terms.
    if(length(linterms > 0)){
      for(l in 1:length(linterms)){
        Xtest[,linterms[l]] = (Xtest[,linterms[l]]-mul[l])/sdl[l]
      }
    }




    ##combine to data frame
    if(linear==FALSE){
      if(intercept==TRUE) {
        X_test = cbind(rep(1, times = dim(Xrt)[1]), Xrt)
      }else{X_test = Xrt}
    } else {
      if(intercept==TRUE) {
        X_test = cbind(rep(1, times = dim(Xrt)[1]), Xtest[,linterms], Xrt)
      }else{
        X_test = cbind(Xtest[,linterms], Xrt)
      }
    }

    if(ytransform == "log") {
      predDist = apply(beta, 1, function(x)exp(as.matrix(X_test)%*%(x) + muy))
      pred = apply(predDist, 1, mean)
      er = sqrt(mean((pred-ytest)^2))
    } else if(is.numeric(y)) {
      pred = ((as.matrix(X_test)%*%(bhat)) + muy)
      er = sqrt(mean((pred-ytest)^2))
    } else {
      predDist = apply(beta, 1, function(x)1/(1+exp(-(as.matrix(X_test)%*%(x)))))
      pred = apply(predDist, 1, mean)
      er = 1-mean(ifelse(pred<0.5, 0, 1)==(as.numeric(ytest)-1))
    }
    out = list(postmean=bhat,postdraws=list(beta=hsmodel[[1]], sigma=hsmodel[[2]], tau=hsmodel[[3]]), rules=rulesFin, df=Xt, y=y, prior=prior, modelstuff=list(linterms=linterms, sdr=sdr, mur=mur, sdl=sdl, mul=mul, ytransform=ytransform, muy=muy,sdy=sdy, intercept=intercept, inputtype=inputtype, model.formula = model.formula, terms_test = terms_test), pred=pred, error=er, X=X)
    class(out) = "HorseRulemodel"
  }
  out
}
