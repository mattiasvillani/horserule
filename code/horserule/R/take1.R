
take1 = function(len) {
  out = c()
  i = 0
  while (i < len){
    out = c(out, i+sample(1:2))
    i = i+2
  }
  out = out[1:len]
  out[seq(1, len, 2)]
}
