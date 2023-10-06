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
(files = list.files('~/Projects/03_MICOM/data/raw/martin8/', full.names = T))


###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: ADAMTS7
# cell: SMC
###########################################################


infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin8//8martin_A7-M1_proteins.csv"
head(read.csv(infile))
bait = 'ADAMTS7'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin8')

# chcek if bait in data
det <- detect(infile, bait); any(det)
read.csv(infile)[det, ]

if (do.write){
  data <- prepare(bait = 'A7', infile = infile,
                  filter.ignore = 'ADAMTS7', verbose = T, raw = T,
                  cols = c('Accession', 'Intensity.TMT6.126', 'Intensity.TMT6.128',
                           'Intensity.TMT6.127', 'Intensity.TMT6.129'))$data
  #data[detect(data, '9000000019'),]$Accession = 'ADAMTS7-3xFLAG'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '8', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}

###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: ADAMTS7
# cell: EC
###########################################################

infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin8//8martin_HAEC-A7_proteins.csv"
head(read.csv(infile))
bait = 'ADAMTS7'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin8')

# chcek if bait in data
det <- detect(infile, bait); any(det)
read.csv(infile)[det, ]

if (do.write){
  data <- prepare(bait = 'A7', infile = infile,
                  filter.ignore = 'ADAMTS7', verbose = T, raw = T,
                  cols = c('Accession', 'Intensity.TMT6.126', 'Intensity.TMT6.130',
                           'Intensity.TMT6.127', 'Intensity.TMT6.131'))$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '8', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: HDAC9
# cell: EC
###########################################################


infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin8//8martin_EC-HDAC_proteins.csv"
head(read.csv(infile))
bait = 'HDAC9'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin8')

# chcek if bait in data
det <- detect(infile, bait); any(det)
read.csv(infile)[det, ]

if (do.write){
  data <- prepare(bait = 'HDAC9', infile = infile,
                  verbose = T, raw = T,
                  cols = c('Accession', 'Intensity.TMT6.128', 'Intensity.TMT6.129',
                           'Intensity.TMT6.130', 'Intensity.TMT6.131'))$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '8', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: KSR2 (I)
# cell: SMC
###########################################################

infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin8//8martin_KSR2_proteins.csv"
head(read.csv(infile))
bait = 'KSR2'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir('KSR2_EC', run = 'martin8')

# chcek if bait in data
det <- detect(infile, bait); any(det)
read.csv(infile)[det, ]

if (do.write){
  data <- prepare(bait = bait, infile = infile,
                  filter.ignore = 'ADAMTS7', verbose = T, raw = T,
                  cols = c('Accession', 'Intensity.TMT6.126', 'Intensity.TMT6.130',
                           'Intensity.TMT6.127', 'Intensity.TMT6.131'))$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '8', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}

