#' RuleHeatmap
#'
#' Produces a heatmap that allows to identify what observations are covered by the most important decision rules. Details can be found in Nalenz & Villani (2017).
#'@import grDevices
#'@import RColorBrewer
#'@param model list containing a model of class "HorseRuleFit".
#'@param k number of most important rules to be shown in the RuleHeat plot.



#' @export

ruleheat = function(model, k){
  Xt = model$df
  Xt[Xt <0] = 0
  Xt[Xt >0] = 1
  postmean = model$postmean
  y = model$y
  imp = abs(postmean)
  imp[model$modelstuff$linterms] = 0
  bhat = postmean[order(imp, decreasing = T)[1:k]]
  Xtemp = Xt[,order(imp, decreasing=T)[1:k]]
  name=c()
  for ( i in 1:k){
    name[i] = paste(i)
  }
  colnames(Xtemp) = name
  bnew= bhat[bhat!=0]
  colcol = c()
  colcol[bnew <0]=9
  colcol[bnew >0]=4
  palette = colorRampPalette(c("white", "aquamarine3"))(2)
  marker = colorRampPalette(brewer.pal(11,"Spectral"))(10)
  if(is.factor(y)==T){
    heatmap(Xtemp, RowSideColors=c("black", "red")[as.numeric(y)], ColSideColors = marker[colcol],
            col = palette, scale='none')
  } else {
    nonz = (y+abs(min(y)))
    vals = nonz/max(nonz)
    int = ((1:11)/11)
    heatmap(as.matrix(Xtemp), RowSideColors=c(grey(rev(1:12)/12)[findInterval(vals, int)+1]), ColSideColors = marker[colcol],
            col = palette, scale='none')
  }
}
