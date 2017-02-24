
createX = function(X, rules, t, corelim=1){
  Xr = matrix(0, nrow=dim(X)[1], ncol=length(rules))
  for (i in 1:length(rules)){
    Xr[eval(parse(text = rules[i])),i] = 1
  }

  Nr = dim(Xr)[2]
  ind = 1:Nr
  if(dim(X)[1]<200){
  t= 0.05
  }
  sup = apply(Xr, 2, mean)
  elim = which((sup<t)|(sup>(1-t)))

  if(length(elim)>0){
    ind = ind[-elim]
  }

  C = cor(Xr[,ind])

  diag(C) = 0
  Nr = dim(Xr[,ind])[2]
  elim=c()
  for(i in 1:(Nr-1)){
    elim = c(elim, which(round(abs(C[i,(i+1):Nr]), digits=4)>=corelim) +i)
  }

  ind = ind[-elim]
  Xr = Xr[,ind]
  rules = rules[ind]
  list(data.matrix(Xr), rules)
}
