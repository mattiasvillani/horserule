rule_length = function(rule) {
  splitted = unlist(strsplit(rule, split = "&"))
  length(splitted)
}
