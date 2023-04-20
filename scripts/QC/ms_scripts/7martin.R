
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
(files = list.files('~/Projects/03_MICOM/data/raw/martin7/', full.names = T))

###########################################################
# date: 09-Jan-2019
# author: Frederik Heymann Lassen
###########################################################

infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin7//7martin_ADAMTS7_proteins.csv"
head(read.csv(infile))
bait = 'ADAMTS7'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin7')

# chcek if bait in data
det <- detect(infile, bait); any(det) #

if (do.write){
  data <- prepare(bait = 'ADAMTS7', cols = c('Accession', 'Intensity.TMT6.126', 'Intensity.TMT6.128','Intensity.TMT6.127', 'Intensity.TMT6.129'),
                  infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '7', bait = bait, cell=cell, bait.vs = 'mock', note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 09-Jan-2019
# author: Frederik Heymann Lassen
# bait: KSR2
# control: control
# cell: SMC
# note: martin claims bait is in the data, but i can't seem to find it
###########################################################

infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin7//7martin_KSR2_proteins.csv"
head(read.csv(infile))
bait = 'KSR2'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin7')

# chcek if bait in data
det <- detect(infile, '3x'); any(det) #
read.csv(infile)[detect(infile, bait), ]

if (do.write){
  data <- prepare(c('SMC', 'KSR'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '7', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}

