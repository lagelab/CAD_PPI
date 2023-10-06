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

(files = list.files('~/Projects/03_MICOM/data/raw/martin9/', full.names = T))

###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: HDAC9
# cell: EC
# note: so  many isoforms. What is going on here. Discuss.
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_HDAC9(N-tag)HAEC_proteins.csv"
head(read.csv(infile))
bait = 'HDAC9' #FN1
cell = 'EC'

# create directories for storaring data
det <- detect(infile, 'HDAC9') # Wow! so many isoforms..
read.csv(infile)[det,]
dirs <- mkdir(bait = bait, cell = cell, run = 'martin9')

if (do.write){
  data <- prepare(bait = c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '9000000018')$data
  #data[detect(data, '9000000018'),]$Accession = '3xFLAG-HDAC9'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}

###########################################################
# date: 10-Jan-2019
# author: Frederik Heymann Lassen
# bait: KIAA1462
# cell: EC
# note: so  many isoforms. What is going on here. Discuss.
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_EC-KIAA1462-KIAA1_proteins.csv"
head(read.csv(infile))
bait = 'KIAA1462' #FN1
cell = 'EC'

# create directories for storaring data
det <- detect(infile, 'FLAG'); any(det)
read.csv(infile)[det,]
dirs <- mkdir(bait = bait, cell = cell, run = 'martin9')

if (do.write){
  data <- prepare(c(cell, 'KIAA'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = 'JCAD', cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 10-Jan-2019
# author: Frederik Heymann Lassen
# bait: KIAA1462
# cell: SMC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_SMC-KIAA1462-KIAA1_proteins.csv"
head(read.csv(infile))
bait = 'KIAA1462'
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'KIAA1462'); any(det)
read.csv(infile)[det,]
dirs <- mkdir(bait = bait, cell = cell, run = 'martin9')

if (do.write){
  data <- prepare(c(cell, 'KIAA'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = 'JCAD', cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}

###########################################################
# date: 10-Jan-2019
# author: Frederik Heymann Lassen
# bait: KSR2
# cell: EC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_EC_KSR2_proteins.csv"
head(read.csv(infile))
bait = 'KSR2'
cell = 'EC'

# create directories for storaring data
det <- detect(infile, 'KSR'); any(det)
read.csv(infile)[det,]
dirs <- mkdir(bait = bait, cell = cell, run = 'martin9')

if (do.write){
  data <- prepare(c(cell, 'KSR'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}

###########################################################
# date: 10-Jan-2019
# author: Frederik Heymann Lassen
# bait: FLT1
# cell: SMC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_SMC_FLT1_proteins.csv"
head(read.csv(infile))
bait = 'FLT1'
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'FLT'); any(det)
read.csv(infile)[det,]
dirs <- mkdir(bait = bait, cell = cell, run = 'martin9')

#
if (do.write){
  data <- prepare(c(cell, 'FLT'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}

###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (WT)
# control: control
# cell: SMC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_ARHRGF(WTmut)-HCASMC_proteins.csv"
head(read.csv(infile))
bait = 'WT' #FARHGEF26
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'ARHGEF')
read.csv(infile)[det,]
dirs <- mkdir(name='ARHGEF26(WT)_SMC', run = 'martin9')

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '3xFLAG-ARHGEF26')$data
  #data[detect(data,'900'),]$Accession = '3xFLAG-ARHGEF26'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = 'ARHGEF26', cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}

###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (MT)
# control: control
# cell: SMC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_ARHRGF(WTmut)-HCASMC_proteins.csv"
head(read.csv(infile))
bait = 'VL29L' #FARHGEF26
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'ARHGEF26')
read.csv(infile)[det,]
dirs <- mkdir(name='ARHGEF26(MT)_SMC', run = 'martin9')

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '3xFLAG-ARHGEF26')$data
  #data[detect(data,'900'),]$Accession = '3xFLAG-ARHGEF26'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = 'ARHGEF26(MT)', cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 08-Jan-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (MT)
# control: control
# cell: SMC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_ARHRGF(WTmut)-HCASMC_proteins.csv"
head(read.csv(infile))
bait = 'WT' #FARHGEF26
cell = 'SMC'

# create directories for storaring data
det <- detect(infile, 'ARHGEF26')
read.csv(infile)[det,]
dirs <- mkdir(name='ARHGEF26(WTandMT)_SMC', run = 'martin9')

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '3xFLAG-ARHGEF26', control = 'VL29L')$data
  #data[detect(data,'900'),]$Accession = '3xFLAG-ARHGEF26'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = 'ARHGEF26', cell=cell, bait.vs = 'ARHGEF26(MT)', note = '3xFlag')), quote = F, sep = '\t')
}



############################################################
############################################################
############################################################

###########################################################
# date: 09-Jan-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (WT)
# control: control
# cell: teloHAEC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_ARHRGF(STmut)-teloHAEC_proteins.csv"
head(read.csv(infile))
bait = 'WT' #FARHGEF26
cell = 'Telo'

# create directories for storaring data
det <- detect(infile, 'ARHGEF26')
read.csv(infile)[det,]
dirs <- mkdir(name='ARHGEF26(WT)_Telo', run = 'martin9')

#
if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '9000000002|9000000001')$data
  #data[detect(data,'002'),]$Accession = '3xFLAG-ARHGEF26-MT'
  #data[detect(data,'001'),]$Accession = '3xFLAG-ARHGEF26-WT'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = 'ARHGEF26', cell='teloHAEC', bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 09-Jan-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (MT)
# control: control
# cell: teloHAEC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_ARHRGF(STmut)-teloHAEC_proteins.csv"
head(read.csv(infile))
bait = 'V29L' #FARHGEF26
cell = 'Telo'

# create directories for storaring data
det <- detect(infile, 'ARHGEF26')
read.csv(infile)[det,]
dirs <- mkdir(name='ARHGEF26(MT)_Telo', run = 'martin9')

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = '9000000002|9000000001')$data
  #data[detect(data,'002'),]$Accession = '3xFLAG-ARHGEF26-MT'
  #data[detect(data,'001'),]$Accession = '3xFLAG-ARHGEF26-WT'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = 'ARHGEF26(MT)', cell='teloHAEC', bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}

###########################################################
# date: 09-Jan-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (WT)
# control: ARHGEF26 (MT)
# cell: teloHAEC
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_ARHRGF(STmut)-teloHAEC_proteins.csv"
head(read.csv(infile))
bait = 'WT' #FARHGEF26
cell = 'Telo'
dirs <- mkdir(name='ARHGEF26(WTandMT)_Telo', run = 'martin9')

if (do.write){
  data <- prepare(c(cell, bait), control = 'V29L', infile = infile, verbose = T, raw = T, filter.ignore = '9000000002|9000000001')$data
  #data[detect(data,'002'),]$Accession = '3xFLAG-ARHGEF26-MT'
  #data[detect(data,'001'),]$Accession = '3xFLAG-ARHGEF26-WT'
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = 'ARHGEF26', cell='teloHAEC', bait.vs = 'ARHGEF26(MT)', note = '3xFlag')), quote = F, sep = '\t')
}



############################################################
############################################################
############################################################



###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (MT)
# control: control
# cell: EC
# note: Bait not present despite being listed in master
# list as present.
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_ARHRGF(WTmut)-HAEC_proteins.csv"
head(read.csv(infile))
bait = 'WT' #FARHGEF26
cell = 'EC'

# create directories for storaring data
det <- detect(infile, 'flag|FLAG|ARHGEF26|SGEF|HMFN1864|CSGEF')
read.csv(infile)[det,]

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = 'FLAG',
                  cols = c('Accession', 'Intensity.EC_WT1.TMT6.126.', 'Intensity.EC_Mock1.TMT6.130.',
                           'Intensity.EC_WT2.TMT6.127.', 'Intensity.EC_Mock2.TMT6.131.'))$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = 'ARHGEF26', cell='EC', bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}




###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (MT)
# control: control
# cell: EC
# note: Bait not present despite being listed in master
# list as present.
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_ARHRGF(WTmut)-HAEC_proteins.csv"
head(read.csv(infile))
bait = 'MT' #FARHGEF26
cell = 'EC'

# create directories for storaring data
det <- detect(infile, 'ARHGEF26|SGEF|HMFN1864|CSGEF')
read.csv(infile)[det,]

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = 'FLAG',
                  cols = c('Accession', 'Intensity.EC_V29L.1.TMT6.128.', 'Intensity.EC_Mock1.TMT6.130.',
                           'Intensity.EC_V29L.2.TMT6.129.', 'Intensity.EC_Mock2.TMT6.131.'))$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = 'ARHGEF26(MT)', cell='EC', bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26 (WT)
# control: ARHGEF26 (MT)
# cell: EC
# note: Bait not present despite being listed in master
###########################################################

# data loading
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin9//9martin_ARHRGF(WTmut)-HAEC_proteins.csv"
head(read.csv(infile))
bait = 'MT' #FARHGEF26
cell = 'EC'

# create directories for storaring data
det <- detect(infile, 'ARHGEF26|SGEF|HMFN1864|CSGEF')
read.csv(infile)[det,]

if (do.write){
  data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, filter.ignore = 'FLAG',
                  cols = c('Accession', 'Intensity.EC_WT1.TMT6.126.', 'Intensity.EC_V29L.1.TMT6.128.',
                           'Intensity.EC_WT2.TMT6.127.', 'Intensity.EC_V29L.2.TMT6.129.'))$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #®colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '9', bait = 'ARHGEF26', cell='EC', bait.vs = 'ARHGEF26(MT)', note = '3xFlag')), quote = F, sep = '\t')
}



