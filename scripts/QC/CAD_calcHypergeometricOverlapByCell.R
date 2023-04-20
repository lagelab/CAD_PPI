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

# write out
#write.table(data[['EC_only']], 'data/interactor_lists/run056/COMBINED_EC_only.InteractorTable.txt', quote = F, row.names = F, sep = '\t')
#write.table(data[['SMC_only']], 'data/interactor_lists/run056/COMBINED_SMC_only.InteractorTable.txt', quote = F, row.names = F, sep = '\t')
#write.table(data[['EC_SMC_only']], 'data/interactor_lists/run056/COMBINED_EC_SMC_only.InteractorTable.txt', quote = F, row.names = F, sep = '\t')

data <- data[21:26]
lapply(data, function(x){
 paste(sum(x$significant),'+',sum(!x$significant))
})

#data <- data[21:26]
reactome <- fread('derived/tables/210422_reactome_msigdb_c2_table.csv')

# databases
database <- list(
  reactome = reactome[,c(3,2)],
  msigdb_h = msigdb_h_table,
  bp = goa_bp_table[,c(1,3)],
  cc = goa_cc_table[,c(1,3)],
  mf = goa_mf_table[,c(1,3)]
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
  go_table <- database[[dbname]]
  colnames(go_table) <- c('genes','geneset')
  terms <- unique(go_table$geneset)
  head(go_table)
  
  # database background
  db_background <- unique(sort(go_table$genes))
  
  result_conditional <- lapply(sort(names(data)), function(sheet){
    
    # get bait and print current data
    bait <- get_bait(sheet)
    print(sheet)
    
    go_enrichment <- do.call(rbind, lapply(terms, function(term){

      # get data (currently GLOBAL analysis)
      df <- as.data.frame(data[[sheet]])
      df_background <- data.frame(gene = db_background[! db_background %in% df$gene], significant = F)
      df_global <- rbind(df[,c('gene','significant')], df_background)
      
      # get term
      go_df <- go_table[go_table$geneset %in% term,]
      go_df$significant <- TRUE
      go_df <- go_df[,c('genes','significant')]
      
      # calculate overlap
      hypergeom <- calc_hyper(df, go_df, bait = bait, intersectDf = data.frame(intersectN = F))
      outdf <- hypergeom$statistics
      outdf$dataset <- term
      outdf$comparison <- 'Conditional Enrichment. FDR < 0.1 & LogFC > 0 (IntersectN = F) '
      outdf$enrichment<- 'conditional'
      outdf$overlap_genes <- paste(hypergeom$genes$mylist$successInSample_genes, collapse = ';')
      
      # calculate global enrichmeent
      hypergeom_global <- calc_hyper(df_global, go_df, bait = bait, intersectDf = data.frame(intersectN = F))
      outdf_global <- hypergeom_global$statistics
      outdf_global$dataset <- term
      outdf_global$comparison <- 'Global Enrichment. FDR < 0.1 & LogFC > 0 (IntersectN = F)'
      outdf_global$enrichment<- 'global'
      outdf_global$overlap_genes <- paste(hypergeom_global$genes$mylist$successInSample_genes, collapse = ';')

      # combine data
      outdf_combined <- rbind(outdf, outdf_global)
      return(outdf_combined)
      
    }))
    
    # summarize data
    go_enrichment$FDR <- NA
    go_enrichment$FDR[go_enrichment$enrichment== 'conditional'] <-  stats::p.adjust(go_enrichment$pvalue[go_enrichment$enrichment== 'conditional'], method = 'fdr')
    go_enrichment$FDR[go_enrichment$enrichment!= 'conditional'] <-  stats::p.adjust(go_enrichment$pvalue[go_enrichment$enrichment!= 'conditional'], method = 'fdr')
    go_enrichment <- go_enrichment[,c(1,8,9,2:7,12,10,11)]
    go_enrichment <- go_enrichment[order(go_enrichment$pvalue), ]
    go_enrichment$list_name <- sheet
    
    return(go_enrichment)
  })
  
  # organize data
  names(result_conditional) <- sort(names(data))
  order <- c(allnames[21:26], allnames[1:20])
  result_conditional <- result_conditional[order]
  result_conditional <- null_omit(result_conditional)
  
  # write file out
  outfile_xlsx = paste0('derived/210422_',dbname,'_network_only_table.xlsx')
  writexl::write_xlsx(result_conditional,outfile_xlsx)

}


#graphics.off()


# read tables and generate heatmaps('derived/tables/)
#pdf('~/Projects/14_micom_clean/MICOM/derived/plots/210409_GO_enrichment_terms.pdf', width = 18, height = 15)
#files <- list.files('derived/tables/', pattern = '210409',  full.names = T)[-1]
#files <- files[!grepl('global',files)]
#for (f in files){
# in_sheets <- readxl::excel_sheets(f)
# d <- lapply(in_sheets,function(x) as.data.frame(readxl::read_xlsx(f, sheet = x)))
# names(d) <- in_sheets
# title = paste(basename(f),'- Geneset Enrichment Analysis')
# ylab = paste0('Gene Ontology (',unlist(lapply(strsplit(f, split = '_'), index, 2)), ')')
# plt <- make_cell_heatmap(d, main = title, ylab = ylab)
# print(plt)
#}
#graphics.off()












