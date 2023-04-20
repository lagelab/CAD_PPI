
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
(files = list.files('~/Projects/03_MICOM/data/raw/martin5/', full.names = T))

###########################################################
# date: 15-dec-2019
# author: Frederik Heymann Lassen
# bait: BCAS3
# cell: SMC
# bait recovered: yes
# note: manually renamed bait row.
###########################################################

infile = files[4]
head(read.csv(infile))
bait = 'BCAS'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(name='BCAS3_SMC', run = 'martin5')

# chcek if bait in data
detect(infile, bait, T)
det <- detect(infile, bait); any(det)
read.csv(infile)[det, ] # symbol is in data but not following convention

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '9000000009')$data
  #data[detect(data, '9000000009'),]$Accession = '3xFLAG-BCAS3'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '5', bait = 'BCAS3', cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}




###########################################################
# date: 15-dec-2019
# author: Frederik Heymann Lassen
# bait: KCNK5
# cell: SMC
# bait recovered: yes
###########################################################

infile = files[8]
head(read.csv(infile))
bait = 'KCNK'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(name='KCNK5_SMC', run = 'martin5')

# chcek if bait in data
detect(infile, bait, T) # symbol is in data but not following convention

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '9000000010')$data
  #data[detect(data, '9000000010'),]$Accession = 'KCNK5-3xFLAG'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '5', bait = 'KCNK5', cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: FLT1
# cell: EC
# bait recovered: bo
###########################################################

infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin5//5martin_FLT_proteins.csv"
head(read.csv(infile))
bait = 'FLT1'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(name='FLT1_EC', run = 'martin5')

# chcek if bait in data
detect(infile, bait, T)
det <- detect(infile, bait); any(det)

if (do.write){
  data <- prepare(c(cell, 'FLT'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '5', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}

###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: KIAA1462
# cell: EC
# bait recovered: Yes (under alias: JCAD)
###########################################################

infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin5//5martin_KIAA1462_proteins.csv"
head(read.csv(infile))
bait = 'KIAA1462'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(name='KIAA1462(JCAD)_EC', run = 'martin5')

# chcek if bait in data
detect(infile, 'KIAA1462|JCAD', T) # found as JCAD, but only one unique protein
det <- detect(infile, bait); any(det)

if (do.write){
  data <- prepare(c(cell, 'KIAA'), infile = infile, verbose = T, raw = T, filter.ignore = 'JCAD')$data
  #data[detect(data, 'JCAD'),]$Accession = 'KIAA1462'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '5', bait = 'JCAD', cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 13-dec-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26
# cell: THP1
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin5//5martin_ARHGEF28_proteins.csv"
head(read.csv(infile))
bait = 'ARHGEF26(WT)'
cell = 'THP1'

# create directories for storaring data
dirs <- mkdir(name = 'ARHGEF26_THP1', run = 'martin5')

# Check data for bait
det <- detect(infile, '^ARHGEF26$'); any(det)

if (do.write){
  data <- prepare(c('THP', 'A28'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '5', bait = 'ARHGEF26', cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}


