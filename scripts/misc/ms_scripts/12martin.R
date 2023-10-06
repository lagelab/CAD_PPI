## initial analysis
if (F){
  .libPaths('~/Toolbox/rlib')
  setwd('~/Projects/03_MICOM/')
  library(pRroteomics)
  library(dplyr)
}
(files = list.files('~/Projects/03_MICOM/data/raw/martin12/', full.names = T))

###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: EDN1 (recovered)
# cell: EC
# Note: FLAG protein needs to be renamed..
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin12//12martin_EDN1_NoNorm_First_UniprotIsoTrembl_proteins.csv"
head(read.csv(infile))
bait = 'EDN1'
cell = 'EC'

# create directories for storaring data
det <- detect(infile, 'EDN1') # bait sucessfully found (1 unique peptide)
read.csv(infile)[det,]
dirs <- mkdir(bait = bait, cell = cell, run = 'martin12')

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = bait)$data
  #data[detect(data,'9000000025'),]$Accession = 'EDN1-3xFLAG'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '12', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}


###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: FLT1 (not recovered)
# cell: SMC
# note: failed IP.
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin12//12martin_FLT1_NoNorm_First_UniprotIsoTrembl_proteins.csv"
head(read.csv(infile))
bait = 'FLT1'
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'FLT1') # not recovered
read.csv(infile)[det,]
dirs <- mkdir(bait = bait, cell = cell, run = 'martin12')

if (do.write){
  data <- prepare(c(cell, 'FLT'), infile = infile, verbose = T, raw = T, filter.ignore = bait)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '12', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}

###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: JCAD (not recovered)
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin12//12martin_JCAD_NoNorm_First_UniprotIsoTrembl_proteins.csv"
head(read.csv(infile))
bait = 'JCAD'
cell = 'EC'

# create directories for storaring data
det <- detect(infile, 'JCAD|KIAA1462') # not recovered
read.csv(infile)[det,]
dirs <- mkdir(bait = bait, cell = cell, run = 'martin12')

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = bait)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '12', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}

###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: PHACTR1 (recovered)
# cell: SMC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin12//12martin_PHACTR1_NoNorm_First_UniprotIsoTrembl_proteins.csv"
head(read.csv(infile))
bait = 'PHACTR1' # PHACTR1
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'PHACTR1')
read.csv(infile)[det,]
dirs <- mkdir(bait = bait, cell = cell, run = 'martin12')

if (do.write){
  data <- prepare(c(cell, 'PHACTR'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '12', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
}

