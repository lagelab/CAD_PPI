setwd('~/Projects/03_MICOM/')
devtools::load_all('~/Projects/08_genoppi_lagelab/Genoppi/')

paths <- list.files('~/Projects/03_MICOM/Downloads/23JUL20_interactor_list/', full.names = T)

goa <- list('goa_bp' = goa_bp_table, 'goa_cc' = goa_cc_table, 'goa_mf' = goa_mf_table)

msig <- list(
  'msigdb_c1' = msigdb_c1_table,                
  'msigdb_c2' = msigdb_c2_table,                
  'msigdb_c3' = msigdb_c3_table,                
  'msigdb_c4' = msigdb_c4_table,                
  'msigdb_c5' = msigdb_c5_table,                
  'msigdb_c6' = msigdb_c6_table,                
  'msigdb_c7' = msigdb_c7_table,                
  'msigdb_ch' = msigdb_h_table,
  'hgnc' = hgnc_group_table
)


pdf('~/Desktop/24JUL20_geneset_enrichment.pdf', width = 20, height = 20)
for (fname in paths){
  
  print(fname)
  df <- read.csv(fname, sep = '\t')
  df <- df[,c('GeneName', 'Interactor')]
  colnames(df) <- c('gene', 'significant')
  
  for (pname in names(goa)){
    
    print(pname)
    ref <- goa[[pname]]
    colnames(ref) <- c('gene','id','group')
    ref$significant <- T
    
    enrichment <- suppressWarnings(calc_adjusted_enrichment(df, ref, col.by = 'group', intersectN = F))
    #print(head(enrichment[,c('list_name','pvalue','BH.FDR')]))
    
    outdf <- enrichment
    outdf$successInSampleGenes <- NULL
    outdf$successInSampleGenes <- enrichment$successInSampleGenes
    outdf <- outdf[outdf$BH.FDR <= 0.25,]
    write.table(outdf, paste0('~/Desktop/',pname,'_',basename(fname)), sep = '\t', quote = F, row.names = F)
    
    ## plot data
    outdf$log10pvalue <- -log10(outdf$pvalue)
    outdf$log10fdr <- -log10(outdf$BH.FDR)
    
    p2 <- plot_tissue_enrichment(outdf, 
                                 'list_name', 
                                 'log10fdr', 
                                 pvalue.line = -log10(0.05),
                                 xlab = paste('Geneset:',pname),
                                 ylab = '-log10(FDR)'
    ) + ggtitle(strsplit(basename(fname), split = '\\.')[[1]][1])
    print(p2)
    
  }
  
  for (pname in names(msig)){
    
    print(pname)
    ref <- msig[[pname]]
    colnames(ref) <- c('gene','group')
    ref$significant <- T
    
    enrichment <- suppressWarnings(calc_adjusted_enrichment(df, ref, col.by = 'group', intersectN = F))
    #print(head(enrichment[,c('list_name','pvalue','BH.FDR')]))
    
    outdf <- enrichment
    outdf$successInSampleGenes <- NULL
    outdf$successInSampleGenes <- enrichment$successInSampleGenes
    outdf <- outdf[outdf$BH.FDR <= 0.25,]
    write.table(outdf, paste0('~/Desktop/',pname,'_',basename(fname)), sep = '\t', quote = F, row.names = F)
    
    ## plot data
    outdf$log10pvalue <- -log10(outdf$pvalue)
    outdf$log10fdr <- -log10(outdf$BH.FDR)
    
    p2 <- plot_tissue_enrichment(outdf, 
                                 'list_name', 
                                 'log10fdr', 
                                 pvalue.line = -log10(0.05),
                                 xlab = paste('Geneset:',pname),
                                 ylab = '-log10(FDR)'
    ) + ggtitle(strsplit(basename(fname), split = '\\.')[[1]][1])
    print(p2)
    
  }
}
graphics.off()


# collect data and save in excel file
paths = list.files('~/Desktop/micom_geneset_analysis/', full.names = T)
groups = c('^(?!.*SMC).*EC','^(?!.*EC).*SMC','EC-SMC') # perl regex
group_names = c('EC','SMC','EC-SMC')
dfnames = c("list_name","successInSample_count", "sample_count", "notSample_count",   
    "success_count","population_count", "pvalue", "BH.FDR", "successInSampleGenes")

lapply(1:3, function (i){
  
  # hack the data.frames to contain whitespace
  group = groups[i]
  grouped_paths <- paths[grepl(group,paths, perl = T)]
  lst <- lapply(grouped_paths, read.csv, sep = '\t', stringsAsFactors = F)
  lst <- lapply(lst, function(x) as.data.frame(rbind(x, rep('', ncol(x)))))
  lst <- lapply(lst, function(x){
    colnames(x) <- dfnames
    return(x)})
  
  names(lst) <- unlist(lapply(strsplit(basename(grouped_paths), '_COMBINED'),function(x) x[1]))
  mat <- as.data.frame(do.call(rbind, lst))
  outname = group_names[i]
  write.table(mat, paste0('~/Desktop/24JUL20_micom_geneset_',outname,'_network_analysis.tsv'), sep = '\t', quote = F, row.names = T)
  
})


