
if (F){
  ## initial analysis
  .libPaths('~/Toolbox/rlib')
  setwd('~/Projects/03_MICOM/')
  library(pRroteomics)
  library(dplyr)
}

(files = list.files('~/Projects/03_MICOM/data/raw/martin15/', full.names = T))

###########################################################
# date: 15-dec-2019
# author: Frederik Heymann Lassen
# bait: EDNRA1
# cell: EC
# Note: bait could not be recovered. No interactors.
###########################################################

# data loading
infile = files[2]
head(read.csv(infile))
bait = 'EDNRA'
cell = 'EC'

# create directories for storaring data
any(detect(infile, 'EDNRA')) # <-- bait not detected
dirs <- mkdir(bait = bait, cell = cell, run = 'martin15')
data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T)$data
write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '15', bait = bait, cell=cell, bait.vs = 'Mock', note = 'endogenous')), quote = F, sep = '\t')

###########################################################
# date: 15-dec-2019
# author: Frederik Heymann Lassen
# bait: PLPP3
# cell: EC
###########################################################

# data loading
infile = files[5]
head(read.csv(infile))
bait = 'PLPP3'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin15')

data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T)$data
write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '15', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
#data <- data[,c('Accession','rep1','rep2')]
#data <- data[complete.cases(data), ]



###########################################################
# date: 15-dec-2019
# author: Frederik Heymann Lassen
# bait: FLT1
# cell: EC
###########################################################
set.seed(1)

# data loading
infile = files[3]
head(read.csv(infile))
bait = 'FLT'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin15')
data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T)$data
write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '15', bait = 'FLT1', cell=cell, bait.vs = 'Mock', note = 'endogenous')), quote = F, sep = '\t')



###########################################################
# date: 15-dec-2019
# author: Frederik Heymann Lassen
# bait: FLT1
# cell: SMC
###########################################################
set.seed(1)

# data loading
infile =  "/Users/flassen/Projects/03_MICOM/data/raw/martin15//15martin_FLT1-HCASMC_NoNorm_First_proteins.csv"
head(read.csv(infile))
bait = 'FLT'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin15')
det <- detect(infile, 'FLT')
read.csv(infile)[det, ]


data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = 'FLT1')$data
write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '15', bait = bait, cell=cell, bait.vs = 'Mock', note = 'endogenous')), quote = F, sep = '\t')
#data <- data[,c('Accession','rep1','rep2')]
#data <- data[complete.cases(data), ]



###########################################################
# date: 16-dec-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (WT)
# cell: SMC
# experiment TMP
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin15//15martin_ARHGEF(WT&mut)_NoNorm_First_proteins.csv"
head(read.csv(infile))
bait = 'WT'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(name = 'ARHGEF26(WT)_SMC', run = 'martin15')

# Check data for bait
det <- detect(infile, 'ARHGEF26'); any(det) # true, bait is there.
read.csv(infile)[det,] # looks like only one unique was found..

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '9000000002')$data
  #data[detect(data, '9000000002'),]$Accession = '3xFLAG-ARHGEF26'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '15', bait = 'ARHGEF26', cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
  #data <- data[,c('Accession','rep1','rep2')]
  #data <- data[complete.cases(data), ]
}


###########################################################
# date: 16-dec-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (MT)
# cell: SMC
# experiment TMP
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin15//15martin_ARHGEF(WT&mut)_NoNorm_First_proteins.csv"
head(read.csv(infile))
bait = 'V29'
cell = 'SMC'

x <- read.csv(infile)
pairs(log10(x[,6:9]))
# create directories for storaring data
dirs <- mkdir(name = 'ARHGEF26(MT)_SMC', run = 'martin15')

# Check data for bait
det <- detect(infile, 'ARHGEF26|FLAG'); any(det) # true, bait is there.
read.csv(infile)[det,] # looks like only one unique was found..

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '9000000002')$data
  #data[detect(data, '9000000002'),]$Accession = '3xFLAG-ARHGEF26-MT'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '15', bait = 'ARHGEF26(MT)', cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
}


###########################################################
# date: 16-dec-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (WT\MT)
# cell: SMC
# control: wildtype
# experiment TMP
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin15//15martin_ARHGEF(WT&mut)_NoNorm_First_proteins.csv"
head(read.csv(infile))
bait = 'WT'
cell = 'SMC'

# create directories for storaring data
dirs <- mkdir(name = 'ARHGEF26(WTandMT)_SMC', run = 'martin15')

# Check data for bait
det <- detect(infile, 'ARHGEF26'); any(det) # true, bait is there.
read.csv(infile)[det,] # looks like only one unique peptide was found is MS..

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '9000000002',
                  cols = c('Accession', 'Intensity.SMC.WT.1.TMT6.126.', 'Intensity.SMC.V29L.1.TMT6.128.',
                           'Intensity.SMC.WT.2.TMT6.127.', 'Intensity.SMC.V29L.2.TMT6.129.'))$data
  #data[detect(data, '9000000002'),]$Accession = '3xFLAG-ARHGEF26'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '15', bait = 'ARHGEF26', cell=cell, bait.vs = 'ARHGEF26(MT)', note = '3xFLAG')), quote = F, sep = '\t')
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
}





