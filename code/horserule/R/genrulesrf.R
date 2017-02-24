#' @importFrom randomForest randomForest
#' @importFrom randomForest combine
#' @import inTrees

genrulesRF = function(X, y, nt,S ,L){
  N = dim(X)[1]
  sf = min(1, (11*sqrt(N)+1)/N)
  mn = 2+floor(rexp(1, 1/(L-2)))
  ns = S
  forest = randomForest(x = X, y=y, sampsize = sf*N ,replace=F, ntree =1, maxnodes=mn, nodesize = ns)
  for(i in 2:nt) {
    mn = 2+floor(rexp(1, 1/(L-2)))
    ns = S
    model1 = randomForest(x = X, y=y, sampsize = sf*N ,replace=F, ntree =1, maxnodes=mn, nodesize = ns)
    forest = combine(forest, model1)
  }
  treelist = RF2List(forest)
  rules = extractRules(treeList=treelist, X=X, ntree=nt, maxdepth=15)
  rules = c(rules)
  rules = rules[take1(length(rules))]
  rulesmat = matrix(rules)
  colnames(rulesmat) = "condition"
  metric = getRuleMetric(rulesmat,X,y)
  pruned = pruneRule(metric, X, y, 0.025, typeDecay = 1)
  unique(pruned[,4])
}
