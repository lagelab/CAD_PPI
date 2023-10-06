
if (F){
  .libPaths('~/Toolbox/rlib')
  setwd('~/Projects/03_MICOM/')
  library(pRroteomics)
  library(dplyr)
}
(files = list.files('~/Projects/03_MICOM/data/raw/martin6/', full.names = T))

###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: HDAC9
# control: control
# cell: SMC
###########################################################

infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin6//6martin_HDAC9_proteins.csv"
head(read.csv(infile))
bait = 'HDAC9'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin6')

# chcek if bait in data
det <- detect(infile, 'flag'); any(det)
read.csv(infile)[det, ]

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '6', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}



###########################################################
# date: 13-dec-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26(WT)
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin6//6martin_ARHGEF26(WT-mutant)_proteins.csv"
head(read.csv(infile))
bait = 'ARHGEF26(WT)'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(name = 'ARHGEF26(WT)_EC', run = 'martin6')

# Check data for bait
det <- detect(infile, 'FLAG'); any(det)
read.csv(infile)[det,]

if (do.write){
  data <- prepare(c('EC', 'WT'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '6', bait = 'ARHGEF26', cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}

###########################################################
# date: 13-dec-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26(MT)
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin6//6martin_ARHGEF26(WT-mutant)_proteins.csv"
head(read.csv(infile))
bait = 'ARHGEF26(MT)'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(name = 'ARHGEF26(MT)_EC', run = 'martin6')

# Check data for bait
det <- detect(infile, '^ARHGEF26$'); any(det) # true, bait is there.
read.csv(infile)[det,] # looks like only one unique was found..

if (do.write){
  data <- prepare(c('EC', 'V29L'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '6', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}



###########################################################
# date: 16-dec-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (WTandMT)
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin6//6martin_ARHGEF26(WT-mutant)_proteins.csv"
head(read.csv(infile))
bait = 'ARHGEF26(MT)'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(name = 'ARHGEF26(WTandMT)_EC', run = 'martin6')

# Check data for bait
det <- detect(infile, '^ARHGEF26$'); any(det) # true, bait is there.
read.csv(infile)[det,] # looks like only one unique was found..

if (do.write){
  data <- prepare(c('EC', 'WT'), control = 'V29L', infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '6', bait = 'ARHGEF26', cell=cell, bait.vs = 'ARHGEF26(MT)', note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 13-dec-2019
# author: Frederik Heymann Lassen
# bait: ADAMTS7
# cell: SMC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin6//6martin_ADAMTS7_proteins.csv"
head(read.csv(infile))
bait = 'ADAMTS7'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(name = 'ADAMTS7_SMC', run = 'martin6')

# Check data for bait
det <- detect(infile, 'ADAM'); any(det)

if (do.write){
  data <- prepare(bait, infile = infile, cols = c('Accession', 'Intensity.iTRAQ4.114', 'Intensity.iTRAQ4.116',
                                                  'Intensity.iTRAQ4.115', 'Intensity.iTRAQ4.117'), verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '6', bait = bait, cell=cell, bait.vs = 'mock', note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 13-dec-2019
# author: Frederik Heymann Lassen
# bait: FLT1
# cell: SMC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin6//6martin_FLT1_proteins.csv"
head(read.csv(infile))
bait = 'FLT1'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(name = 'FLT1_SMC', run = 'martin6')

# Check data for bait
det <- detect(infile, 'FLT'); any(det)

if (do.write){
  data <- prepare(c('SMC','FLT'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '6', bait = bait, cell=cell, bait.vs = 'mock', note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 13-dec-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26
# cell: THP1
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin6//6martin_ARHGEF26_proteins.csv"
head(read.csv(infile))
bait = 'ARHGEF26'
cell = 'THP1'

# create directories for storaring data
dirs <- mkdir(name = 'ARHGEF26_THP', run = 'martin6')

# Check data for bait
det <- detect(infile, 'FLT'); any(det)

if (do.write){
  data <- prepare(c('THP','A26'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '6', bait = bait, cell=cell, bait.vs = 'mock', note = '3xFlag')), quote = F, sep = '\t')
}


