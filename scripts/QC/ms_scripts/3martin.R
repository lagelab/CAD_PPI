if (F){
  ## initial analysis
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
(files = list.files('~/Projects/03_MICOM/data/raw/martin3/', full.names = T))

###########################################################
# date: 19-dec-2019
# author: Frederik Heymann Lassen
# bait: BCAS3
# cell: EC
# bait recovered: YES
###########################################################


infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin3//3martin_BCAS3_proteins.csv"
head(read.csv(infile))
bait = 'BCAS3'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin3')

# chcek if bait in data
det <- detect(infile, bait); any(det)
read.csv(infile)[det, ]

# Write
if (do.write){
  data <- prepare(bait = c(cell,bait), infile = infile, filter.ignore = c('0004'), verbose = T, raw = T)$data
  #data[grepl('0004',data$Accession),]$Accession = '3xFLAG-BCAS3'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '3', bait = bait, bait.vs = 'Mock', cell=cell, note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 19-dec-2019
# author: Frederik Heymann Lassen
# bait: KCNK5
# cell: EC
# bait recovered: YES
###########################################################


infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin3//3martin_KCNK5_proteins.csv"
head(read.csv(infile))
bait = 'KCNK5'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin3')

# chcek if bait in data
det <- detect(infile, bait); any(det)
read.csv(infile)[det, ]

# Write
if (do.write){
  data <- prepare(bait = c(cell,bait), infile = infile, filter.ignore = c('9000000005'), verbose = T, raw = T)$data
  #data$Accession[1] = c('KCNK5-3xFLAG')
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip('Whitehead', martin = '3', bait = bait, bait.vs = 'Mock', cell=cell, note = '3xFlag')), quote = F, sep = '\t')
}

