read_omim <- function(regex1 = 'coronary', regex2='susceptibility', path = '~/Projects/15_genesets/genesets/data/omim/genemap2.txt'){
  d <- suppressWarnings(fread(path)) # the end of the file contains a readme, which reuslts in a warning..
  rows <- apply(d, 1, function(x) tolower(paste(x, collapse = ' ')))
  bool <- grepl(tolower(regex1), rows) & grepl(tolower(regex2), rows)
  d <- d[bool,]
  if (nrow(d) > 0){
    d$gene <- d$`Approved Symbol`
    d$gene[d$gene == ''] <- unlist(lapply(strsplit(d$`Gene Symbols`[d$gene == ''], split = ', '), function(x) x[1]))
    d$significant = TRUE
    d <- d[!duplicated(d$gene),]
    return(d)
  } else {
    return(NULL)
  }
}
