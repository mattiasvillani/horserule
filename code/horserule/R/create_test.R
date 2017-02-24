
createXtest = function(X, rules) {
  Xr = matrix(0, nrow=dim(X)[1], ncol=length(rules))
  for (i in 1:length(rules)){
    Xr[eval(parse(text = rules[i])),i] = 1
  }
  data.matrix(Xr)
}
