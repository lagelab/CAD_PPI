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
(files = list.files('~/Projects/03_MICOM/data/raw/martin11/', full.names = T))

###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: FLT1 (not recovered)
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin11//11martin_FLT1_NoNorm_First_UniprotIsoTrembl_proteins.csv"
head(read.csv(infile))
bait = 'FLT1' #FN1
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'FLT')
dirs <- mkdir(bait = bait, cell = cell, run = 'martin11')

#
if (do.write){
  data <- prepare(c(cell, 'FLT'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '11', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}


###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: KSR2 (not recovered)
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin11//11martin_KSR2_NoNorm_First_UniprotIsoTrembl_proteins.csv"
head(read.csv(infile))
bait = 'KSR2' #FN1
cell = 'EC'

# create directories for storaring data
det <- detect(infile, 'KSR|FLAG'); any(det)
dirs <- mkdir(bait = bait, cell = cell, run = 'martin11')

#
if (do.write){
  data <- prepare(c(cell, 'KSR'), infile = infile, verbose = T, raw = T,
                  cols = c('Accession', 'Intensity.iTRAQ4.114', 'Intensity.iTRAQ4.116',
                           'Intensity.iTRAQ4.115', 'Intensity.iTRAQ4.117'))$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '11', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}


###########################################################
# date: 10-Jan-2019
# author: Frederik Heymann Lassen
# bait: FN1 (recovered as alias:FINC)
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin11//11martin_FN1_NoNorm_First_UniprotIsoTrembl_proteins.csv"
head(read.csv(infile))
bait = 'FN1' #FN1
cell = 'EC'

# create directories for storaring data
det <- detect(infile, 'CIG|MSF|FNZ|LETS|GFND|FINC') # Alias is FINC...
read.csv(infile)[det,]
dirs <- mkdir(bait = bait, cell = cell, run = 'martin11')

#
if (do.write){
  data <- prepare(c(cell, 'FN'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '11', bait = bait, cell=cell, bait.vs = 'Mock', note = 'Endogenous')), quote = F, sep = '\t')
}

###########################################################
# date: 10-Jan-2019
# author: Frederik Heymann Lassen
# bait: EDNRA (not recovered)
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin11//11martin_EC-EDNRA1_NoNorm_First_UniprotIsoTrembl_proteins.csv"
head(read.csv(infile))
bait = 'EDNRA'
cell = 'EC'

# create directories for storaring data
det <- detect(infile, 'EDNRA'); any(det)
dirs <- mkdir(bait = bait, cell = cell, run = 'martin11')

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '11', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}


###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: EDNRA (not recovered)
# cell: SMC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin11//11martin_SMC-EDNRA1_NoNorm_First_UniprotIsoTrembl_proteins.csv"
head(read.csv(infile))
bait = 'EDNRA'
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'EDNRA'); any(det)
dirs <- mkdir(bait = bait, cell = cell, run = 'martin11')

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '11', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}

