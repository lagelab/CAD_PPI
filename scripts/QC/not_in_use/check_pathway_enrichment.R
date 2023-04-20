
files <- list.files('~/Projects/03_MICOM/derived/', pattern = '14JAN.*\\.txt$',
                    recursive = T, full.names = T)
for (f in files){
  tab <- read.table(f, header = T)
  tab$pathway <- unlist(lapply(tab$gene, function(x) paste(lst[[x]], collapse = ';' )))
  tab$pathway <- (lapply(tab$gene, function(x) lst[[x]]))
  
  # 
  all_pathways <- unlist(tab[tab$FDR < 0.1 & tab$logFC > 0, ]$pathway)
  print(table(all_pathways))
  
  #
  #fname = gsub('\\.txt','_W_PATHWAYS\\.txt', f)
  #write.table(tab, fname)
  }


kegg <- lst
save(kegg, file = 'kegg_pathways.RData')

files <- list.files('~/Projects/03_MICOM/derived/', pattern = '(14|15)JAN.*\\.txt$',
                    recursive = T, full.names = T)
for (f in files){
  tab <- read.table(f, header = T)
  tab$pathway <- unlist(lapply(tab$gene, function(x) paste(lst[[x]], collapse = ';' )))
  #tab$pathway <- (lapply(tab$gene, function(x) lst[[x]]))
  outname <- gsub('\\.txt', 'W_PATHWAYS.txt', f)
  t <- data.frame(gene = tab$gene, pathway = tab$pathway)
  write.table(t, outname, quote = F, row.names = F, sep = '\t')
}

