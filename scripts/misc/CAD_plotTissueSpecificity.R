# gtex analysis of enrichment
# by frhl
# data: 21-04-13

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
gtex_cat <- fread('~/Projects/15_genesets/genesets/data/gtex/GTEX.tstat.categories.genoppi.csv')
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

# manually import networks
read_network <- function(d){return(data.frame(gene = d$GeneName, significant = d$Interactor))}

EC <- read_network(fread('data/interactor_lists/run056/COMBINED_EC.InteractorTable.txt'))
SMC <- read_network(fread('data/interactor_lists/run056/COMBINED_SMC.InteractorTable.txt'))
EC_SMC <- read_network(fread('data/interactor_lists/run056/COMBINED_EC-SMC.InteractorTable.txt'))

# EC/SMC networks
data[['EC']] = EC
data[['SMC']] = SMC
data[['EC_SMC_union']] = EC_SMC

# Tissue Specific networks
data[['EC_only']] = subset(EC, (significant & !(gene %in% SMC$gene[SMC$significant])) | !significant)
data[['SMC_only']] = subset(SMC, (significant & !(gene %in% EC$gene[EC$significant])) | !significant)
data[['EC_SMC_intersect']] = subset(EC_SMC, (significant & (gene %in% EC$gene[EC$significant]) & (gene %in% SMC$gene[SMC$significant])) | !significant)

# subset to only networks
data <- data[-1:-20]

# what data query?
database <- list(
  gtex_rna = gtex_rna
  #gtex_protein = gtex_protein,
  #hpa_rna = hpa_rna
)


allnames <- names(data)

# setup helpers
get_bait <- function(x){unlist(strsplit(x, split = '_'))[1]}


#outfile_pdf = paste0('derived/plots/210409_enrichment_heatmaps.pdf')
#pdf(outfile_pdf, width = 16, height = 13)
## actual script for getting enrichments
for (dbname in names(database)){
  
  print(paste('####', dbname, Sys.time(),'####'))

  # load
  tissue_table <- database[[dbname]]
  terms <- as.character(unique(tissue_table$tissue))
  bonf <- 0.05 / length(terms)
  print(paste(dbname,'- 0.05 / ',  length(terms), '- bonf:', bonf))
  head(tissue_table)
  
  # database background
  db_background <- as.character(unique(sort(tissue_table$gene)))
  
  result_conditional <- lapply(sort(names(data)), function(sheet){
    
    # get bait and print current data
    bait <- get_bait(sheet)
    print(sheet)
    
    tissue_enrichment <- do.call(rbind, lapply(terms, function(term){
      
      # get data (currently GLOBAL analysis)
      df <- as.data.frame(data[[sheet]])
      df_background <- data.frame(gene = db_background[! db_background %in% df$gene], significant = F)
      df_global <- rbind(df[,c('gene','significant')], df_background)
      
      # get term
      tissue_df <- tissue_table[tissue_table$tissue %in% term,]
      
      # calculate overlap
      hypergeom <- calc_hyper(df, tissue_df, bait = bait, intersectDf = data.frame(intersectN = F))
      outdf <- hypergeom$statistics
      #outdf$significant <- outdf$pvalue < bonf
      outdf$dataset <- term
      outdf$comparison <- 'Conditional Enrichment. FDR < 0.1 & LogFC > 0 (IntersectN = F) '
      outdf$enrichment <- 'conditional'
      #outdf$overlap_genes <- paste(hypergeom$genes$mylist$successInSample_genes, collapse = ';')
      
      # calculate global enrichmeent
      hypergeom_global <- calc_hyper(df_global, tissue_df, bait = bait, intersectDf = data.frame(intersectN = F))
      outdf_global <- hypergeom_global$statistics
      #outdf_global$significant <- outdf_global$pvalue < bonf
      outdf_global$dataset <- term
      outdf_global$comparison <- 'Global Enrichment. FDR < 0.1 & LogFC > 0 (IntersectN = F)'
      outdf_global$enrichment <- 'global'
      #outdf_global$overlap_genes <- paste(hypergeom_global$genes$mylist$successInSample_genes, collapse = ';')
      
      # combine data
      outdf_combined <- rbind(outdf, outdf_global)
      
      return(outdf_combined)
      
    }))
    
    # summarize the data
    tissue_enrichment$FDR <- NA
    tissue_enrichment$FDR[tissue_enrichment$enrichment == 'conditional'] <-  stats::p.adjust(tissue_enrichment$pvalue[tissue_enrichment$enrichment == 'conditional'], method = 'fdr')
    tissue_enrichment$FDR[tissue_enrichment$enrichment != 'conditional'] <-  stats::p.adjust(tissue_enrichment$pvalue[tissue_enrichment$enrichment != 'conditional'], method = 'fdr')
    tissue_enrichment <- tissue_enrichment[,c(1,8,9,2:7,10,11)]
    tissue_enrichment$list_name <- sheet
    #tissue_enrichment <- tissue_enrichment[tissue_enrichment$enrichment == 'conditional',]

    # add tissue categories
    tissue_enrichment <- merge(gtex_cat, tissue_enrichment, by.x = 'Tissue.genoppi', by.y = 'dataset')
    tissue_enrichment$Tissue.genoppi <- NULL
    tissue_enrichment$Number.of.samples <- NULL
    
    # order by p-value
    tissue_enrichment$significant <- tissue_enrichment$pvalue < bonf
    tissue_enrichment <- tissue_enrichment[order(tissue_enrichment$pvalue), ]
    
    return(tissue_enrichment)
  })
  
  # organize data
  names(result_conditional) <- sort(names(data))
  #order <- c(allnames[21:26], allnames[1:20])
  #result_conditional <- result_conditional[order]
  result_conditional <- null_omit(result_conditional)

  
  # write file to out table
  outfile_xlsx = paste0('derived/210422_',dbname,'_all_network_table.xlsx')
  writexl::write_xlsx(result_conditional,outfile_xlsx)
  
}




