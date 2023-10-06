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
(files = list.files('~/Projects/03_MICOM/data/raw/martin10/', full.names = T))

###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: FN1 (recovered as alias:FINC)
# cell: SMC
# note: three different isoforms are detected, but they all have the exact same intensity
# values.
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin10//10martin_SMC_FN1_proteins.csv"
head(read.csv(infile))
bait = 'FN' #FN1
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'CIG|MSF|FNZ|LETS|GFND|FINC') # Alias is FINC...
read.csv(infile)[det,]  # three isoforms.. But they all have some intensity values.
dirs <- mkdir(bait = bait, cell = cell, run = 'martin10')

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '!FINC')$data
  #data[detect(data,'!FINC'),]$Accession = 'FN1'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '10', bait = 'FN1', cell=cell, bait.vs = 'Mock', note = 'endogenous')), quote = F, sep = '\t')
}

###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: JCAD (recovered)
# cell: SMC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin10//10martin_SMC_JCAD_proteins.csv"
head(read.csv(infile))
bait = 'JCAD'
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'JCAD|KIAA1462')
read.csv(infile)[det,]
dirs <- mkdir(bait = bait, cell = cell, run = 'martin10')

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '9000000024')$data
  #data[detect(data, '9000000024'),]$Accession = 'JCAD-3xFLAG'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '10', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}


###########################################################
# date: 10-Jan-2019
# author: Frederik Heymann Lassen
# bait: EDNRA (not recovered)
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin10//10martin_ED_EDNRA1_proteins.csv"
head(read.csv(infile))
bait = 'EDNRA'
cell = 'EC'

# create directories for storaring data
det <- detect(infile, 'EDNRA'); any(det)
dirs <- mkdir(bait = bait, cell = cell, run = 'martin10')

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '10', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}


###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: EDNRA (not recovered)
# cell: SMC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin10//10martin_SMC_EDNRA1_proteins.csv"
head(read.csv(infile))
bait = 'EDNRA'
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'EDNRA'); any(det)
dirs <- mkdir(bait = bait, cell = cell, run = 'martin10')

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '10', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}


###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: FLT1 (not recovered)
# cell: SMC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin10//10martin_SMC_FLT1_proteins.csv"
head(read.csv(infile))
bait = 'FLT1'
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'FLT'); any(det)
dirs <- mkdir(bait = bait, cell = cell, run = 'martin10')

if (do.write){
  data <- prepare(c(cell, 'FLT'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '10', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}

