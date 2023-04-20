description.keyword <- function(keyword = 'GN=', from = df$Description){
  
  x1 = strsplit(from, split = keyword)
  x2 = unlist(lapply(x1, function(x) x[length(x)]))
  x3 = strsplit(x2, split = '\\ ')
  x4 = unlist(lapply(x3, function(x) x[1]))
  return(x4)
}

