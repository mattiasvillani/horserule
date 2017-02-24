calc_prior = function(rules, Xr,alpha, beta) {
  sup = apply(Xr, 2, function(x)mean(x>0))
  prior = c()
  len = unlist(lapply(rules, rule_length))
  sup = sqrt((sup)*(1-sup))
  for(i in 1:length(rules)){
    prior[i] = ((2*sup[i])^(alpha))/(len[i]^beta)
  }
  prior
}
