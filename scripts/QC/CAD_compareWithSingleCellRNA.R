setwd('~/Projects/03_MICOM/')
library(ggrepel)
library(gplots)
library(ggplot2)
library(RColorBrewer)
devtools::load_all('~/Projects/08_genoppi_lagelab/Genoppi/')


# ready paths and prefixes
date = '201120'
run = 'run056'

# setup paths and extract tier 1 paths
prefix = paste0(date,'_',run,'_micom')
files = read.csv('run056_filtered_tier1_paths.tsv', sep ='\t')$x
baits = unlist(lapply(strsplit(files, split = '\\.|vs'), function(x)x[3]))
cells = unlist(lapply(strsplit(files, split = '\\.|vs'), function(x)x[5]))

# organize data into matrix by bait x cell
mylist = lapply(unique(baits), function(bait){
  
  bool.ec = cells == 'EC' & baits == bait
  bool.smc = cells == "SMC" & baits == bait
  ec.files = files[bool.ec][1] # for now, remove duplicates
  smc.files = files[bool.smc][1] # for now, remove duplicates
  if (!any(bool.ec)) ec.files <- NA
  if (!any(bool.smc)) smc.files <- NA
  return(data.frame(bait = bait, EC = ec.files, SMC = smc.files))
  
})

# now we have contrast in each cell type
contrast = do.call(rbind, mylist)
contrast = contrast[!as.logical(apply(contrast, 1, function(x) sum(is.na(x)))),] # remove NA's for now


###################
# prep scRNA data #
###################

SE <- function(x) sd(x)/sqrt(length(x))
tstat <- function(x) { (mean(x) - x)/SE(x)}

if (F){
  
  scRNA <- read.csv('data/external/rna_single_cell_type 2.tsv', sep = '\t') #human protein atlas RNA seq
  scRNA.cells <- unique(scRNA$Cell.type)
  scRNA.genes <- unique(scRNA$Gene.name)
  
  # for each gene, compute a tstat for specific expression in cortex versus all non-brain samples
  
  gene = 'KCNK5'
  cell = 'Mucus-secreting cells'
  
  
  cell.tstat <- tstat(scRNA$NX[scRNA$Cell.type == cell])
  
  cell.tstat[scRNA$Gene.name[scRNA$Cell.type == cell] == gene]
  
  #scRNA.sig <-lapply(scRNA.cells, function(cell){
  #  print(cell)
  #  subdf <- scRNA[scRNA$Cell.type == cell,]
  #  subdf$significant <- subdf$NX >= quantile(subdf$NX, 0.9)
  #  return(subdf)
  #})
  #scRNA <- do.call(rbind, scRNA.sig)
  #colnames(scRNA) <- c('ensemble_gene_id', 'gene', 'cell', 'NX', 'significant')
 
  
}




#################################
# for each pair, get some stats #
# ###############################


scRNA <- read.table('data/external/hpa_rna_single_cell_type_sig.txt', sep = '\t')
pdf('test_enrichment_pvalue.pdf', width = 8, height = 7)
for (i in 1:nrow(contrast)){
  
  # setup paths
  bait = contrast$bait[i]
  pathEC <- contrast$EC[i]
  pathSMC <- contrast$SMC[i]
  nameEC = tools::file_path_sans_ext(basename(pathEC))
  nameSMC = tools::file_path_sans_ext(basename(pathSMC))
  
  # setup data
  dfEC = read.csv(pathEC, row.names = NULL, sep = '\t')
  dfSMC = read.csv(pathSMC, row.names = NULL, sep = '\t')
    
  # hypergeometric
  resEC <- calc_adjusted_enrichment(dfEC, scRNA, col.by = 'cell', verbose = T)
  resSMC <- calc_adjusted_enrichment(dfSMC, scRNA, col.by = 'cell', verbose = T)

  # combine
  cond1 = 'EC'
  cond2 = 'SMC'
  pos = resEC[,c(1,7)]
  pos$condition = as.factor('EC')
  neg = resSMC[,c(1,7)]
  neg$condition = as.factor('SMC')
  
  # setup data.frame
  dt <- rbind(pos, neg)
  colnames(dt) <- c('list_name', 'FDR', 'condition')
  dt$cell <- factor(dt$list_name, levels = unique(dt$list_name[order(dt$FDR)]))
  dt$list_name <- NULL
  dt$FDR <- -log10(dt$FDR)
  
  xmi <- -10
  xma <- 10
  
  plt = ggplot(data = dt, aes(x = reorder(cell, -FDR), fill = condition)) +
    geom_bar(stat = "identity", data = subset(dt, condition == cond1), aes(y=FDR)) +
    geom_bar(stat = "identity", data = subset(dt, condition == cond2), aes(y=FDR * (-1)) ) +
    scale_y_continuous(limits = c(xmi, xma), breaks = seq(xmi, xma, 10), labels = abs(seq(xmi, xma, 10))) + 
    theme(axis.text = element_text(colour = "black")) + 
    coord_flip() + 
    ggtitle(paste('scRNA data'))+
    ylab("-log10(Hypergeometric P-value)") + 
    xlab("") +
    theme_minimal()
  print(plt)
  
}
graphics.off()






