
whitehead.accession.uniprot.map <- function(x, table = uniprot_to_gene){
  ids = unlist(strsplit(as.character(x), split = '\\|'))
  ids = unique(unlist(lapply(strsplit(ids, split = '\\-'), function(y) y[1])))
  ids_mapped = unlist(lapply(ids, function(id) map_uniprot_to_gene_nopaste(id, table)[1]))
  ids_mapped = unique(ids_mapped)
  if (all(is.na(ids_mapped))) return(NA)
  ids_mapped = ids_mapped[!is.na(ids_mapped)]
  genes = paste(ids_mapped, collapse = '|')
  if (genes == '') browser()
  return(genes)
}


#missing_data = lapply(dfs1, function(x) {
#  #genes = x$gene.first.description
#  x = x[x$species == 'HUMAN',]
#  x = x[x$Uncharacterized == F,]
#  x = x[!is.na(x$Accession),]
#  genes = as.character(x$gene.first.description)
#  x[is.na(genes),]$Accession
#})

#unlist(missing_data)


#
#missingness = unlist(lapply(dfs1, function(x) {
#  x = x[x$species == 'HUMAN',]
#  x = x[x$Uncharacterized == F,]
#  x = x[!is.na(x$Accession),]
#  genes = as.character(x$gene.first.description)
#  return(sum(is.na(genes))/length(genes))
#  }))
