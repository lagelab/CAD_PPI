


# gtex ?

# misigdb ?

calc_nested_enrichment <- function(data, genesets, intersectN = F, background_genes = NULL){
  
  if (! all(c('gene','significant','geneset') %in% colnames(genesets))) stop('missing collumns for genesets list')
  
  
  result <- lapply(sort(names(data)), function(sheet){
    
    # get bait and print current data
    browser()
    df <- as.data.frame(data[[sheet]])
    bait <- get_bait(sheet)
    print(sheet)
    
    # Allow to compare with all human genes background
    if (!is.null(background_genes)){
      df_background <- data.frame(gene = df$gene[! background_genes %in% df$gene], significant = F)
      df <- rbind(df[,c('gene','significant')], df_background)
    }
    
    # get terms
    terms <- unique(genesets$geneset)
    go_enrichment <- do.call(rbind, lapply(terms, function(term){
      
      # get term
      genesets_df <- genesets[genesets$geneset %in% term,]
      #genesets_df$significant <- TRUE
      genesets_df <- genesets_df[,c('genes','significant')]
      
      # calculate overlap
      hypergeom <- calc_hyper(df, go_df, bait = bait, intersectDf = data.frame(intersectN = intersectN))
      outdf <- hypergeom$statistics
      outdf$dataset <- term
      outdf$comparison <- 'FDR < 0.1 & LogFC > 0 (IntersectN = F)'
      outdf$overlap_genes <- paste(hypergeom$genes$mylist$successInSample_genes, collapse = ';')
      return(outdf)
      
    }))
    
    # summarize data
    go_enrichment$FDR <- stats::p.adjust(go_enrichment$pvalue, method = 'fdr')
    go_enrichment <- go_enrichment[,c(1,8,9,2:7,11,10)]
    go_enrichment <- go_enrichment[order(go_enrichment$pvalue), ]
    go_enrichment$list_name <- sheet
    return(go_enrichment)
  })
  
  
}



calc_nested_enrichment(data, database)

