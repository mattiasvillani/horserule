#' Variable Importance plot
#'
#' Creates a input variable importance plot
#'@param model list containing a model of class "hs_rulefit".
#'@param top If a integer number is given only shows the top most important variables in the plot
#'@param var_names optional vector with the variable names to be shown in plot.
#'@import ggplot2
#'@export

Variable_importance = function(model, top=NULL, var_names=NULL){
  if(class(model)!="HorseRulemodel"){
    stop("Model must be of class HorseRulemodel")
  }

  beta = (model$postdraws)$beta
  rules = model$rules
  inp = dim(model$X)[2]
  lin = (model$modelstuff)$linterms
  start = ifelse(model$modelstuff$intercept ==T, 1, 0) + length(lin)
  p = dim(beta)[2]
  if(is.null(var_names)){
    var_names = sprintf("variable %i", 1:p)
  }
  samp = dim(beta)[1]
  out = matrix(0, nrow=samp, ncol=inp)
  for(i in ((start+1):p)) {
    splitted = unlist(strsplit(rules[i-start], split = " & "))
    len = length(splitted)
    vars = unique(unlist(strsplit(splitted, "in|<=|>"))[seq(from=1, to=len*2, by=2)])
    m = length(vars)
    for(j in 1:m){
      str = vars[j]
      ind = as.numeric(regmatches(str, gregexpr("[0-9]+", str)))
      out[, ind]  = out[, ind] + abs(beta[,i])/m
    }
  }
  if(length(lin)>0){
    for(j in 1:length(lin)){
      out[,lin[j]] = out[,lin[j]] + abs(beta[,j])
    }
  }
  if(!is.null(top)){

    normout = t(apply(out, 1, function(x)((x-min(x))/(max(x)-min(x)))))
    meanimp = apply(normout,2,mean)
    ind = order(meanimp, decreasing = T)[1:top]
    topout = normout[,ind]
    name = c()
    val = c()
    for(j in 1:top){
      name = c(name, rep(var_names[ind][j], times = samp))
      val = c(val, topout[,j])
    }
  } else {
    normout = t(apply(out, 1, function(x)((x-min(x))/(max(x)-min(x)))))
    meanimp = apply(normout,2,mean)
    ind = order(meanimp, decreasing = T)
    topout = normout[,ind]
    name = c()
    val = c()
    for(j in 1:inp){
      name = c(name, rep(var_names[ind][j], times = samp))
      val = c(val, topout[,j])
    }
  }
  frame = data.frame(factor(name, levels=unique(name)), val)
  colnames(frame) =c("Input", "Importance")
  p = ggplot(frame, aes("Input", "Importance"))
  p + stat_boxplot(geom = "errorbar",colour = I("#3366FF"))+
    geom_boxplot(colour = I("#3366FF"))    +theme_bw()  +theme(text = element_text(size=20, angle=90))
}
