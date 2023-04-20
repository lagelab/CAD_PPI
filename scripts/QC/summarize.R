mlist <- lapply(list.files('derived/',recursive = T, pattern = '\\_SUMMARY', full.names = T), function(x) {
  
  res = as.vector(read.csv(x)) 
  if ('X' %in% names(res)) res$X <- NULL
  if (!'data.bait.enrichment' %in% names(res)) res$data.bait.enrichment <- NA

  
  return(res)
  
})



