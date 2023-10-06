
if (F){
  ## initial analysis
  .libPaths('~/Toolbox/rlib')
  setwd('~/Projects/03_MICOM/')
  #source('~/Toolbox/packages/pRoteomics/R/prepare.R')
  devtools::load_all('~/Toolbox/packages/pRoteomics/')
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
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin16//16martin_FLT1_proteins.csv"
bait = 'FLT1'; cell = 'EC'

# create directories for storaring data
any(detect(infile, bait))
read.csv(infile)[detect(infile, bait),]
dirs <- mkdir(bait = bait, cell = cell, run = 'martin16')
data <- prepare(c(cell, 'FLT'), infile = infile, verbose = T, raw = T)$data
write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '16', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')


