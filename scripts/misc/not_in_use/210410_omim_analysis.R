# what genesets are enriched in the WT versus MT.
setwd('~/Projects/14_micom_clean/MICOM/')
library(ggplot2)
library(genoppi)
library(cowplot)
library(writexl)
library(data.table)
source('R/make_cell_heatmap.R')
source('R/read_omim.R')
run = 'run056'

# helpers
index <- function(x, i) x[i]
background_genes <- unique(c(inweb_table$Gene1, inweb_table$Gene2))

# setup paths and run data
files <- list.files('MICOM/data/ms_data_genoppi_input/run056/', full.names = T)
experiments <- read.csv('derived/micom_run056_experiments.csv')
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

# manually import networks
read_network <- function(d){return(data.frame(gene = d$GeneName, significant = d$Interactor))}

EC <- read_network(fread('data/interactor_lists/run056/COMBINED_EC.InteractorTable.txt'))
SMC <- read_network(fread('data/interactor_lists/run056/COMBINED_SMC.InteractorTable.txt'))
EC_SMC <- read_network(fread('data/interactor_lists/run056/COMBINED_EC-SMC.InteractorTable.txt'))

# EC/SMC networks
data[['EC']] = EC
data[['SMC']] = SMC
data[['EC_SMC']] = EC_SMC

# Tissue Specific networks
data[['EC_only']] =  EC[! EC$gene %in% SMC$gene,]
data[['SMC_only']] =  SMC[! SMC$gene %in% EC$gene,]
data[['EC_SMC_intersect']] = EC_SMC[EC_SMC$gene %in% EC$gene & EC_SMC$gene %in% SMC$gene,]





path <- '~/Projects/15_genesets/genesets/data/omim/genemap2.txt'
#re2 <- 'suscep'
#re2 <- 'Dominant'
re2 <- ""

omim <- list(
  stroke = read_omim('stroke', re2, path),
  myocardial = read_omim('myocardial', re2, path),
  hypertension = read_omim('hypertension', re2, path),
  diabetes = read_omim('diabetes', re2, path),
  cardiomyopathy = read_omim('Cardiomyopathy', re2, path),
  vascular = read_omim('Vascular', re2, path)
)

writexl::write_xlsx(omim, 'derived/tables/210410_omim_genesets.xlsx')


#omim <- as.data.frame(do.call(rbind, omim))
#omim <- omim[!duplicated(omim$gene),]
#omim <- list(x = omim)

unlist(lapply(omim, length)) 





## geneset enrichment
result_omim <- lapply(sort(names(data)), function(sheet){
  
  # get bait and print current data
  bait <- get_bait(sheet)
  print(sheet)
  
  omim_enrichment <- do.call(rbind, lapply(names(omim), function(term){
    
    
    # get all data
    df <- as.data.frame(data[[sheet]])
    #df_background <- data.frame(gene = df$gene[! background_genes %in% df$gene], significant = F)
    #df <- rbind(df[,c('gene','significant')], df_background)
    gene_df <- as.data.frame(omim[[term]])
    
    if (nrow(gene_df > 0)){
      
      # calculate overlap
      hypergeom <- calc_hyper(df, gene_df, bait = bait, intersectDf = data.frame(intersectN = F))
      outdf <- hypergeom$statistics
      outdf$dataset <- term
      outdf$comparison <- 'FDR < 0.1 & LogFC > 0 (IntersectN = F)'
      outdf$overlap_genes <- paste(hypergeom$genes$mylist$successInSample_genes, collapse = ';')
      return(outdf)
      
    }
    
    
  }))
  
  # summarize data
  omim_enrichment$FDR <- stats::p.adjust(omim_enrichment$pvalue, method = 'fdr')
  omim_enrichment <- omim_enrichment[,c(1,8,9,2:7,11,10)]
  omim_enrichment <- omim_enrichment[order(omim_enrichment$pvalue), ]
  omim_enrichment$list_name <- sheet
  return(omim_enrichment)
})

result1 <- do.call(rbind, result_omim )
result1 <- result1[order(result1$pvalue),]
result1$FDR <- stats::p.adjust(result1$pvalue)
write_xlsx(result1, 'derived/tables/210410_omim_enrichment.xlsx')
