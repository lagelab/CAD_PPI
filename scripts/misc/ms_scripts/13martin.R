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

(files = list.files('~/Projects/03_MICOM/data/raw/martin13/', full.names = T))

###########################################################
# date: 15-dec-2019
# author: Frederik Heymann Lassen
# bait: PLPP3
# cell: SMC
# bait recovered: yes
# note: There were no interactors for PLPP3 in inweb. Searching
# under different alises did not yield any interactors in inweb.
###########################################################

infile = files[2]
head(read.csv(infile))
bait = 'PLPP3'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin13')

# chcek if bait in data
det <- detect(infile, 'PLPP3'); any(det)
read.csv(infile)[det, ]

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = bait)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '13', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}




###########################################################
# date: 15-dec-2019
# author: Frederik Heymann Lassen
# bait: PHACTR1
# cell: EC
# bait recovered: yes
# note: acession ID for bait (PHACTR1) did not indicate a
# species. For subsequently analysis, I assumed that it was
# a human protein and therefore included in the analysis.
###########################################################

infile = files[1]
head(read.csv(infile))
bait = 'PHACTR1'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin13')

# chcek if bait in data
detect(infile, 'FLAG', T); any(det)
read.csv(infile)[det, ] # acession ID can not be recognized, also no species declaration: gi|9000000026|!PHACTR1-FLAG

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '9000000026')$data
  #data[detect(data, '9000000026'), ]$Accession = 'PHACTR1-FLAG'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '13', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}



