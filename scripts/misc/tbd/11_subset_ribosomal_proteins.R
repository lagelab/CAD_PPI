

# dirs
from <- '/Users/flassen/Projects/03_MICOM/data/genoppi_input/run019'
to <- '/Users/flassen/Projects/03_MICOM/data/genoppi_input/run019_noRP'
qc.ips <- readLines('~/Projects/03_MICOM/21JUL20_micom_qced_ip_paths.txt')

#
from.files <- list.files(from, full.names = T)[list.files(from) %in% qc.ips]
genes.removed <- list()
genes.pct <- list()
genes.removed.sig <- list()
genes.pct.sig <- list()

for (path in from.files){
  
  df <- read.csv(path, sep = '\t')
  gene.RPL <- grepl('^RPL',df$gene)
  gene.RPS <- grepl('^RPS',df$gene)
  gene.keep <- !(gene.RPL | gene.RPS)
  newdf <- df[gene.keep, ]
  outpath <- file.path(to, basename(path))
  #write.csv(newdf, outpath, quote = F, row.names = F)
  print(paste(nrow(df), '->', nrow(newdf), '\t diff:', nrow(df)-nrow(newdf)))
  genes.removed[[path]] <- df$gene[!gene.keep]
  genes.removed.sig[[path]] <- df$gene[!gene.keep & df$significant]
  genes.pct[[path]] <- length(df$gene[!gene.keep])/nrow(df)
  genes.pct.sig[[path]] <- length(df$gene[!gene.keep & df$significant])/nrow(df[df$significant,])
  
}

#write.csv(stack(genes.removed), '~/Desktop/28JUL20_micom_ribosomal_proteins_by_experiment.csv', quote = F)

# what ips are too high?
highlight <- genes.pct.sig[unlist(genes.pct.sig) > 0.20]
highlight <- data.frame(id = basename(names(highlight)), value = round(unlist(highlight), 2))
row.names(highlight) <- NULL
combi <- apply(highlight, 1, paste, collapse = ': ')

# write
pdf('~/Desktop/28JUL20_run019_ribosomal_proteins.pdf', width = 10, height = 10)
par(mfrow=c(2,2))
hist(unlist(lapply(genes.removed, length)), main = 'Absolute count \n (All detected proteins)', xlab = 'ribosomal proteins (n)', xlim = c(0,100), breaks = 10, labels = F)
hist(unlist(genes.pct),  main = 'Relative amount \n (All detected proteins)', xlab = '% ribosomal proteins', xlim = c(0,1), labels = F)
hist(unlist(lapply(genes.removed.sig, length)), main = 'Absolute count \n (FDR <= 0.1)', xlab = 'ribosomal proteins (n)', xlim = c(0,100), breaks = 10, labels = F)
hist(unlist(genes.pct.sig),  main = 'Relative amount \n (FDR <= 0.1)', xlab = '% ribosomal proteins', xlim = c(0,1), labels = F)
text(paste(combi,  collapse = '\n'), x = 0.5, y = 8, cex = 0.4)
graphics.off()

