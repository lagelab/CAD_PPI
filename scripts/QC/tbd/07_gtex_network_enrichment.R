devtools::load_all('~/Projects/08_genoppi_lagelab/Genoppi/')


paths <- list.files('~/Projects/03_MICOM/Downloads/23JUL20_interactor_list/', full.names = T)


pdf('~/Desktop/23JUL20_gtex_enrichment.pdf')
for (fname in paths){
  
  print(fname)
  df <- read.csv(fname, sep = '\t')
  df <- df[,c('GeneName', 'Interactor')]
  colnames(df) <- c('gene', 'significant')
  enrichment <- suppressWarnings(calc_adjusted_enrichment(df, gtex_table))
  print(enrichment[,c('list_name','pvalue','BH.FDR')])

  ## write data
  outdf <- enrichment
  outdf$successInSampleGenes <- NULL
  outdf$successInSampleGenes <- enrichment$successInSampleGenes
  write.table(outdf, paste0('~/Desktop/',basename(fname)), sep = '\t', quote = F, row.names = F)
  
  ## plot data
  outdf$log10pvalue <- -log10(outdf$pvalue)
  outdf$log10fdr <- -log10(outdf$BH.FDR)
  
  
  p1 <- plot_tissue_enrichment(outdf, 
                         'list_name', 
                         'log10pvalue', 
                         pvalue.line = -log10(0.05),
                         xlab = 'Tissue (GTEx)',
                         ylab = '-log10(p-value)'
  ) + ggtitle(strsplit(basename(fname), split = '\\.')[[1]][1])
  
  
  
  p2 <- plot_tissue_enrichment(outdf, 
                         'list_name', 
                         'log10fdr', 
                         pvalue.line = -log10(0.05),
                         xlab = 'Tissue (GTEx)',
                         ylab = '-log10(FDR)'
  ) + ggtitle(strsplit(basename(fname), split = '\\.')[[1]][1])
  print(p1)
  print(p2)
  
}
graphics.off()



