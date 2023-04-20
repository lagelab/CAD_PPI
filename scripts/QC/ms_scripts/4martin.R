
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
(files = list.files('~/Projects/03_MICOM/data/raw/martin4/', full.names = T))

###########################################################
# date: 16-dec-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (WT)
# cell: THP
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin4//4martin_ARHGEF26_proteins.csv"
head(read.csv(infile))
bait = 'ARHGEF26'
cell = 'THP1'

# create directories for storaring data
dirs <- mkdir(name = 'ARHGEF26(WT)_THP', run = 'martin4')

# Check data for bait
det <- detect(infile, 'FLAG'); any(det)
read.csv(infile)[det,]

# Write
if (do.write){
  data <- prepare(c('THP', 'A26'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip('Whitehead', martin = '4', bait = 'ARHGEF26', bait.vs = 'Mock', cell=cell, note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 16-dec-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (WT)
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin4//4martin_ARHGEF26(WT-mutant)_proteins.csv"
head(read.csv(infile))
bait = 'ARHGEF26(WT)'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(name = 'ARHGEF26(WT)_EC', run = 'martin4')

# Check data for bait
det <- detect(infile, '^ARHGEF26$|flag'); any(det)
read.csv(infile)[det,]

# Write
if (do.write){
  data <- prepare(c('EC', 'WT'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip('Whitehead', martin = '4', bait = 'ARHGEF26', bait.vs = 'Mock', cell=cell, note = '3xFlag')), quote = F, sep = '\t')
}

###########################################################
# date: 16-dec-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (MT)
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin4//4martin_ARHGEF26(WT-mutant)_proteins.csv"
head(read.csv(infile))
bait = 'ARHGEF26(MT)'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(name = 'ARHGEF26(MT)_EC', run = 'martin4')

# Check data for bait
det <- detect(infile, '^ARHGEF26$'); any(det) # true, bait is there.
read.csv(infile)[det,] # looks like only one unique was found..

# Write
if (do.write){
  data <- prepare(c('EC', 'mt'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '4', bait = bait, bait.vs = 'Mock', cell=cell, note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 16-dec-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (WTandMT)
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin4//4martin_ARHGEF26(WT-mutant)_proteins.csv"
head(read.csv(infile))
bait = 'ARHGEF26(MT)'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(name = 'ARHGEF26(WTandMT)_EC', run = 'martin4')

# Check data for bait
det <- detect(infile, '^ARHGEF26$'); any(det) # true, bait is there.
read.csv(infile)[det,] # looks like only one unique was found..

# Write
if (do.write){
  data <- prepare(c('EC', 'WT'), control = 'mt', infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip('Whitehead', martin = '4', bait = 'ARHGEF26', bait.vs = 'ARHGEF26(MT)', cell=cell, note = '3xFlag')), quote = F, sep = '\t')
}

###########################################################
# date: 16-dec-2019
# author: Frederik Heymann Lassen
# bait: KSR2
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin4//4martin_KSR2_proteins.csv"
head(read.csv(infile))
bait = 'KSR2'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(name = 'KSR2_EC', run = 'martin4')

# Check data for bait
det <- detect(infile, 'KSR2'); any(det)
read.csv(infile)[det,]

if (do.write){
  data <- prepare(c('EC', 'KSR'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '4', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}

