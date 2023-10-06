
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

(files = list.files('~/Projects/03_MICOM/data/raw/martin14/', full.names = T))


#(files = list.files('~/Projects/03_MICOM/data/genoppi_input/raw/martin14/', full.names = T))

###########################################################
# date: 15-dec-2019
# author: Frederik Heymann Lassen
# bait: EDNRA
# cell: SMC
# bait recovered: yes
# result: failed IP
# description: One unique peptide for the bait can be found
# which has a 100% match with the expected bait (EDNRA) using
# BLAST. Plotting correlations of ratios and intensity values
# seem to indicate that the data is fine. One possibility is that
# there was some kind of mislabelling on the data.
###########################################################

infile = files[2]
head(read.csv(infile))
bait = 'EDNRA'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin14')

# chcek if bait in data
det <- detect(infile, bait); any(det)
read.csv(infile)[det, ] ## only one unique detected in iTRAQ

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = 'EDNRA')$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '14', bait = bait, cell=cell, bait.vs = 'Mock', note = 'endogenous')), quote = F, sep = '\t')
}



###########################################################
# date: 15-dec-2019
# author: Frederik Heymann Lassen
# bait: FLT
# cell: SMC
# bait recovered: no
###########################################################

infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin14//14martin_FLT1_proteins.csv"
head(read.csv(infile))
bait = 'FLT1'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin14')

# chcek if bait in data
det <- detect(infile, 'VGRF1|FLT'); any(det)
read.csv(infile)[det, ]

if (do.write){
  data <- prepare(c(cell, 'FLT'), infile = infile, verbose = T, raw = T, filter.ignore = 'FLT')$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '14', bait = 'FLT1', cell=cell, bait.vs = 'Mock', note = 'endogenous')), quote = F, sep = '\t')
}


###########################################################
# date: 15-dec-2019
# author: Frederik Heymann Lassen
# bait: FLT (FLT or JCAD)
# cell: EC
# bait recovered: no
###########################################################

infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin14//14martin_FLT1-JCAD_proteins.csv"
head(read.csv(infile))
bait = 'FLT1-JCAD'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin14')

# chcek if bait in data
det <- detect(infile, 'FLT|JCAD'); any(det) # JCAD present, FLT1 is not.
read.csv(infile)[det, ]

if (do.write){
  data <- prepare(c(cell, 'FLT'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '14', bait = 'FLT-JCAD', cell=cell, bait.vs = 'Mock', note = 'endogenous')), quote = F, sep = '\t')
}

###########################################################
# date: 15-dec-2019
# author: Frederik Heymann Lassen
# bait: FLT1
# cell: EC
# bait recovered: no
###########################################################

#infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin14//14martin_FLT1_proteins.csv"
#head(read.csv(infile))
#bait = 'FLT1'
#cell = 'SMC'
#
# create directories for storaring data
#dirs <- mkdir(bait = bait, cell = cell, run = 'martin14')
#
# chcek if bait in data
#det <- detect(infile, 'FLT'); any(det)
#read.csv(infile)[det, ]
#
#if (do.write){
#  data <- prepare(c(cell, 'FLT'), infile = infile, verbose = T, raw = T)$data
#  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
#  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
#  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '14', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
#}



