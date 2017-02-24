#' @import gbm
#' @import inTrees


genrulesGBM = function(X, y, nt, S, L) {
  N = dim(X)[1]
  sf = min(1, (11*sqrt(N)+1)/N)
  mn = 2+floor(rexp(1, 1/(L-2)))
  ns = S
  dist = ifelse(is.numeric(y), "gaussian", "bernoulli")
  if (is.numeric(y)==F){
    y = as.numeric(y)-1
  }
  model1 = gbm.fit(x = X, y=y, bag.fraction = sf,n.trees =1, interaction.depth=(mn/2)
                   ,shrinkage = 0.01,distribution = dist, verbose = F, n.minobsinnode = ns)
  for(i in 2:nt) {
    mn = 2+floor(rexp(1, 1/(L-2)))
    model1$interaction.depth = (mn/2)
    model1 = gbm.more(model1, n.new.trees=1, verbose = F)
  }
  treelist = GBM2List(model1, X)
  rules = extractRules(treelist, X=X, ntree=nt, maxdepth=15)
  rules = c(rules)
  rules = rules[take1(length(rules))]
  rulesmat = matrix(rules)
  colnames(rulesmat) = "condition"
  metric = getRuleMetric(rulesmat,X,y)
  pruned = pruneRule(metric, X, y, 0.025, typeDecay = 1)
  unique(pruned[,4])
}
