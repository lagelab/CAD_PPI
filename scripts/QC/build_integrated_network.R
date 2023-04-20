## initial analysis
if (F){
  .libPaths('~/Toolbox/rlib')
  setwd('~/Projects/03_MICOM/')
  library(dplyr)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(hash)
  library(dplyr)
  library(igraph)
  
  #library(rProteomics, lib.loc = '~/Toolbox/rlib/')
  devtools::load_all('~/Toolbox/packages/pRoteomics/')
}

source('scripts/graph_utils.R')

## Load data
files = list.files('data/genoppi_input/run004', full.names = T)
data = lapply(files, function(df) {read.csv(df,sep='\t') %>% designate(FDR<0.1,logFC>0)})
names(data) = basename(files)
data = data[grepl('nonSGS', names(data)) | !grepl('SGS', names(data))]
info = read.csv('10FEB2020_masterlist_analysis_plan.csv')

## subset good IPs

summary_good_ips <- read.table('~/Projects/03_MICOM/tier1_IPs.txt')
#summary = get_micom_ip_summary(data)
#summary_good_ips = as.character(summary[summary$replicate_correlation > 0.6 & summary$bait_found == TRUE & summary$mock == TRUE,]$fname)
good_ips = data[names(data) %in% summary_good_ips$V1]
df = melt_ips(good_ips)

## 
write.table(df, '06MAR2020_micom_tier_1_ips.csv', quote = F)


## Make upset map
outlist = list()
for (cell in c('EC','SMC')){
  for (bait in unique(df$bait)){
    data_bait_cell = df[df$bait %in% bait & df$cell %in% cell & df$significant == TRUE, ]
    if (nrow(data_bait_cell) > 0){
      outlist[[paste0(bait,'  ',cell)]] = data_bait_cell$gene
    }

  }
}
#outlist[['InWeb_InBiomap']] <- hash::keys(inweb_hash)[hash::keys(inweb_hash)%in% df[df$significant == TRUE, ]$gene]
#outlist[['CAD GWAS Risk Genes']] <- read.table('genelist/cad_gwas_genes_prio1_2.txt')[[1]]

## Print to file
library(RColorBrewer)
library(UpSetR)
pdf('plots/06MAR2020_micom_upset_maps.pdf', width = 12, height = 12)
colors = brewer.pal(2,"Reds") #c('red','green','blue','cyan')
plt = upset(fromList(outlist), order.by = "freq", point.size = 3.5, line.size = 1, sets = names(outlist),
            mainbar.y.label = "Overlap", sets.x.label = paste0("Interactors (FDR<0.1, logFC>0)"), keep.order = T,
            sets.bar.color = c(rep(colors[2],8), rep(colors[1],10)))
print(plt)
graphics.off()







