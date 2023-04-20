# mapping uniprot to gene
map_uniprot_to_gene_nopaste <- function(x, table = uniprot_to_gene) unlist(lapply(x, function(i) {
  if (i %in% table$From){
    return(table[as.character(i) == table$From,]$To)
  } else {
    return(NA)
  }
}))

# mapping tablew
#uniprot_to_gene <- read.delim(url("https://www.uniprot.org/mapping/M20200520216DA2B77BFBD2E6699CA9B6D1C41EB28A32AE3.tab"), stringsAsFactors = F)
map_uniprot_to_gene <- function(x, table = uniprot_to_gene) unlist(lapply(x, function(i) {
  if (i %in% table$From){
    return(paste(table[as.character(i) == table$From,]$To, collapse = ' '))
  } else {
    return(NA)
  }
}))