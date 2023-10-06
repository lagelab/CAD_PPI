
if (F){
  ## initial analysis
  .libPaths('~/Toolbox/rlib')
  setwd('~/Projects/03_MICOM/')
  library(pRroteomics)
  library(dplyr)
}

(files = list.files('~/Projects/03_MICOM/data/raw/martin16/', full.names = T))

########################################################### 
# date: 04-Feb-2020
# author: Frederik Heymann Lassen
# bait: FLT1
# cell: EC
# Note: bait could not be recovered. No interactors.
###########################################################

# data loading
#infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin16//16martin_FLT1_proteins.csv"
#head(read.csv(infile))
#bait = 'FLT1'
#cell = 'EC'

# create directories for storaring data
#any(detect(infile, bait)) 
#read.csv(infile)[detect(infile, bait),]
#dirs <- mkdir(bait = bait, cell = cell, run = 'martin16') 

# Prepare data for analysis by getting replicate folds
#df <- prepare(c(cell, 'FLT'), infile = infile, impute = list(shift = -1.8, stdwidth = 0.3), verbose = T)

# get logfoldchange and do moderated ttest
#data <- df %>% mttest()
#write.table(data,file=dirs$txtpath,row.names=F,sep="\t",quote=F)
#pdf(dirs$pdfpath, height=8, width=8)

## plot volcano and scatter plot
#data %>% designate(FDR < 0.1) %>% plotVolcano(bait, title = 'FLT1 vs Control (ES)')
#r = data %>% designate(FDR < 0.1) %>% plotScatter(bait, title = 'FLT1 Replicate Correlation')

# plot inweb interactors
#known.interactors = interactors('FLT1', T)
#data %>% designate(FDR < 0.1) %>% plotOverlap(bait, known.interactors, 'FLT1 vs Control (ES): InWeb Overlap')
#data %>% designate(FDR < 0.1, logFC > 0) %>% plotOverlap(bait, known.interactors, 'FLT1 vs Control (ES): InWeb Overlap')

# plot proxies
#proxies <- expandProxies(data)
#data %>% designate(FDR < 0.1) %>% plotProxies(proxies)
#graphics.off()