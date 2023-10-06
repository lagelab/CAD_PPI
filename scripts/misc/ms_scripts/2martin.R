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

  #library(rProteomics, lib.loc = '~/Toolbox/rlib/')
  devtools::load_all('~/Toolbox/packages/pRoteomics/')
}

(files = list.files('~/Projects/03_MICOM/data/raw/martin2/', full.names = T))

###########################################################
# date: 09-Jan-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin2//2martin_ARHGF26_NoNorm_First_proteins.csv"
head(read.csv(infile))
bait = 'ARHGEF26' #FN1
cell = 'HEK293T'

# create directories for storaring data
det <- detect(infile, '90000')
read.csv(infile)[det,]
#dirs <- mkdir(name='CXCL12beta_EC', run = 'martin2')

# Write
if (do.write){
  data <- prepare(bait = c('WT'), infile = infile, filter.ignore = c('0001','0002'), verbose = T, raw = T)$data
  #data$Accession[1:2] = c('3xFLAG-ARHGEF26-MT','3xFLAG-ARHGEF26-WT')
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip('Whitehead', martin = '2', bait = 'ARHGEF26', bait.vs = 'Mock', cell=cell, note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 09-Feb-2020
# author: Frederik Heymann Lassen
# bait: CXCL12-beta
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin2//2martin_CXCLR_proteins (1).csv"
head(read.csv(infile))
bait = 'CXCL12-beta' #FN1
cell = 'EC'

# create directories for storaring data
det <- detect(infile, 'CX')
read.csv(infile)[det,]
#dirs <- mkdir(name='CXCL12beta_EC', run = 'martin2')

# Write
if (do.write){
  data <- prepare(bait = c('CXCL'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip('Whitehead', martin = '2', bait = bait, bait.vs = 'Mock', cell=cell, note = '3xFlag')), quote = F, sep = '\t')
}

