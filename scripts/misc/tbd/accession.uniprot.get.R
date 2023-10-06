accession.uniprot.get <- function(accession){
  
  # 
  x = accession
  y = accession
  
  xg = grepl('sp',x) # second last items always
  y[xg] = unlist(lapply(strsplit(x[xg], split='\\|'), function(x) x[length(x)-1]))
  
  # these few items that do not conform are NAs
  xg2 = !xg & grepl('gi', x)
  y[xg2] <- NA
  
  # 
  xg3 = !xg2 & grepl('tr', x)
  y[xg3] <- unlist(lapply(strsplit(x[xg3], split='\\|'), function(x) x[2]))
  
  
  #
  xg4 = !xg & !xg2 & !xg3
  y[xg4] <- unlist(lapply(strsplit(x[xg4], split='\\|'), function(x) x[1]))
  
  return(y)
}

uniprot.isoform.get <- function(uniprot){
  
  # isoform seems to mapped wrong
  
  x = strsplit(uniprot, split = '\\-|\\.')
  x = do.call(rbind, x)
  x[,2][x[,1] == x[,2]] <- NA
  x = data.frame(x)
  x$X2 <- as.numeric(x$X2)
  colnames(x) <- c('uniprot','id')
  return(x)
  
}

