# what genesets are enriched in the WT versus MT.
setwd('~/Projects/14_micom_clean/MICOM/')
library(ggplot2)
library(genoppi)
library(cowplot)
library(writexl)
library(data.table)
source('R/make_cell_heatmap.R')
run = 'run056'

# helpers
index <- function(x, i) x[i]
background_genes <- unique(c(inweb_table$Gene1, inweb_table$Gene2))

# setup paths and run data
files <- list.files('MICOM/data/ms_data_genoppi_input/run056/', full.names = T)
experiments <- read.csv('derived/tables/micom_run056_experiments.csv')
experiments$path <- gsub('data/genoppi_input/','data/ms_data_genoppi_input/',experiments$path)
experiments$submission <- lapply(strsplit(basename(experiments$path), split = '\\.'), index, 2)
experiments$id <- paste(experiments$bait,
                        experiments$cell,
                        ifelse(experiments$How == 'Overexpression', 'OE','EN'), 
                        ifelse(experiments$Facility == 'Broad', 'BR','WH'),
                        toupper(experiments$submission), sep = '_')

# load IP-MS data 
data <- lapply(experiments$path, function(x) fread(x))
names(data) <- experiments$id
baits <- unique(experiments$bait)
cells <- unique(experiments$cell)

# Helper functions
get_bait_name <- function(x) return(unlist(strsplit(x, '_'))[1])
get_cell_name <- function(x) return(unlist(strsplit(x, '_'))[2])
get_type <- function(x) return(unlist(strsplit(x, '_'))[3])
get_facility <- function(x) return(unlist(strsplit(x, '_'))[4])

# collapse a list of IPs into a superset
collapse_ips <- function(lst){
  superset = as.data.frame(do.call(rbind, lapply(lst, function(df){df[,c('gene','significant')]})))
  sig_genes = unique(superset$gene[superset$significant])
  superset$both_insig_and_insig = superset$gene %in% sig_genes & superset$significant == FALSE
  superset$significant = superset$significant | superset$both_insig_and_insig
  outdf = superset[,c('gene','significant')]
  outdf = outdf[!duplicated(outdf), ]
  return(outdf)
}

#for (analysis in c('SMC','EC','(SMC)|(EC)')){
for (analysis in c('(SMC)|(EC)')){
  
  experiments_allowed <- experiments[grepl(analysis, experiments$cell),] 
  outlist <- list()
  outlist_overlap <- list()
  
  for (id1 in experiments_allowed$id){
    data1 = data[[id1]]
    outlist[[id1]] = list()
    outlist_overlap[[id1]] = list()
    bait = experiments$bait[experiments$id == id1]
    for(id2 in experiments_allowed$id){
      
      #if (id1 != id2 & id2 == "EDN1_SMC_OE_BR_MB" & id1 == "EDNRA_SMC_EN_WH_M14") browser()
      #if (id1 == id2) browser()
      # get data
      data2 = data[[id2]]
      superset_id1 = collapse_ips(list(data1))
      superset_id2 = collapse_ips(list(data2))
      statistics = calc_hyper(superset_id1, superset_id2, data.frame(intersectN = T), bait = bait)
      statistics$statistics$pvalue_lt_005 <- ifelse(statistics$statistics$pvalue < 0.05, 'Y','N')
      
      #overlap = statistics$statistics$successInSample_count / statistics$statistics$sample_count
      # calculate overlap
      xsig = superset_id1$gene[superset_id1$significant]
      ysig =  superset_id2$gene[superset_id2$significant]
      overlap = intersect(xsig, ysig)
      set1 <- xsig[! xsig %in% overlap]
      set2 <- ysig[! ysig %in% overlap]
      total <- c(overlap, set1, set2)
      
      # what is the intersect of interactors versus the union of 
      overlap_pct <- length(overlap) / length(total)
      outlist_overlap[[id1]][[id2]] = overlap_pct 
      outlist_overlap[[id1]][[id2]] = overlap_pct 
      
      # generate outdf
      comparisondf <- data.frame(bait1 = get_bait_name(id1),cell1 = get_cell_name(id1), type1 = get_type(id1),factility1 = get_facility(id1),
                                 bait2 = get_bait_name(id2),cell2 = get_cell_name(id2), type2 = get_type(id2),factility2 = get_facility(id2))
      overlapdf <- data.frame(overlap = overlap_pct)
      outdf <- cbind(comparisondf, statistics$statistics, overlapdf)
      outlist[[id1]][[id2]] <- outdf
      
      #outlist[[id1]][[id2]] = statistics$statistics # statistics$statistics$pvalue
    }
  }
  
  
  # make table of data
  res <- do.call(rbind, lapply(outlist, function(x) do.call(rbind, x)))
  
  # deal w
  res$duplicated <- duplicated(as.character(do.call(rbind,lapply(strsplit(rownames(res), '\\.'), function(x) paste(sort(x), collapse = '_')))))
  res <- res[!res$duplicated,]
  res$list_name <- NULL
  write.csv(res,'derived/tables/210418_interactor_overlap.csv', quote = F)
  
  
  # Get data and ready for ggplot
  result_pvalue_raw = as.matrix(do.call(rbind, lapply(outlist, unlist)))
  result_pvalue = format(as.matrix(do.call(rbind, lapply(outlist, unlist))),digits=3)
  result_overlap = round(as.matrix(do.call(rbind, lapply(outlist_overlap, unlist))),3)
  
  #result_pvalue_melted = melt(result_pvalue)
  #result_pvalue_melted$FDR <- stats::p.adjust(result_pvalue_melted$value, method = 'fdr')
  #result_overlap_melted = melt(result_overlap, id = c('Var1', 'Var2'))
  #results_combined = cbind(result_pvalue_melted, result_overlap_melted$value)
  #colnames(results_combined) <- c('Var1', 'Var2','value','overlap')
  #results_combined$overlap_pct <- paste0(results_combined$overlap*100, '%')
  
  result_overlap[is.na(result_overlap)] <- 0
  
  

  # make labels
  label_pvalue = paste0('p-value:\n',result_pvalue)
  
  width = 16
  height = 13
  
  if (analysis == '(SMC)|(EC)') {
    analysis <- 'SMC and EC'
    width = 24
    height = 20
    #label_pvalue = ''
  }
  
  
  
  labels = matrix( paste(paste0(result_overlap*100,'%'), label_pvalue, sep="\n"), 
                   nrow=nrow(result_overlap), dimnames=dimnames(result_overlap) )
  
  
  #initiate cols with all black
  cols <- rep('black', nrow(result_overlap))
  cols[grepl('SMC',row.names(result_overlap))] <- 'red'
  rows <- rep('black', ncol(result_overlap))
  rows[grepl('SMC',colnames(result_overlap))] <- 'red'
  
  library(gplots)
  title = paste('Overlap Heatmap', analysis)
  outfile = paste('derived/plots/210415_overlap', analysis, 'heatmap.pdf',collapse = '_')
  pdf(outfile, width = width, height = height)
  palette = colorRampPalette(c("white", "orange"))(n = 100)
  heatmap.2(t(result_overlap),
            dendrogram = 'column',
            cellnote = t(labels),    # same data set for cell labels
            main = title,          # heat map title
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(16,16),     # widens margins around plot
            col=palette,
            colRow = cols,
            colCol = rows
            
            )      
  graphics.off()
  
  
}

