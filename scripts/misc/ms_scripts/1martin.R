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


  write.ip = function(martin, bait, cell, bait.vs = 'Mock', facility = 'Whitehead', note = 'expanded',imputation = 'MinImputed'){
    return(paste0('Genoppi_',facility,'_martin',martin,'_',imputation,'.',bait,'vs',bait.vs,'.',cell,'.',note,'.tsv'))
  }
}
(files = list.files('~/Projects/03_MICOM/data/raw/martin1/', full.names = T))

###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: ARHGEF26
# control: control
# cell: THP1
###########################################################
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin1//Martin_ARHGEF26_proteins.csv"
head(read.csv(infile))
bait = 'ARHGEF26'
cell = 'THP1'

# create directories for storaring data
dirs <- mkdir(bait = bait, cell = cell, run = 'martin1')

# chcek if bait in data
det <- detect(infile, 'flag'); any(det) #

# Write
if (do.write){
  data <- prepare(bait = c('THP', 'A26'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip('Whitehead', martin = '1', bait = 'ARHGEF26', bait.vs = 'Mock', cell=cell, note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: CXCL12-alpha
# control: control
# cell: THP1
###########################################################
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin1//Martin_CXCL12-alpha_proteins.csv"
head(read.csv(infile))
bait = 'CXCL12-alpha'
cell = 'THP1'

# create directories for storaring data
dirs <- mkdir(name = 'CXCL12alpha_THP1', run = 'martin1')

# chcek if bait in data
det <- detect(infile, 'SDF|CXCL12|CXC'); any(det) #

# Write
if (do.write){
  data <- prepare(bait = c('THP', 'CXCL12'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip('Whitehead', martin = '1', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}



###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: CXCL12-beta
# control: control
# cell: EC
###########################################################
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin1//Martin_CXCL12-beta_proteins.csv"
head(read.csv(infile))
bait = 'CXCL12-beta'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(name = 'CXCL12beta_CXCL12', run = 'martin1')

# chcek if bait in data
det <- detect(infile, 'SDF|CXCL12|CXC'); any(det) #

# Write
if (do.write){
  data <- prepare(bait = c('EC', 'CXCL12'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip('Whitehead', martin = '1', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFlag')), quote = F, sep = '\t')
}


###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: PGFD
# control: MOCK
# cell: EC
###########################################################
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin1//Martin_EC-PDGFD_proteins.csv"
head(read.csv(infile))
bait = 'PDGFD'
cell = 'EC'

# create directories for storaring data
dirs <- mkdir(name = 'PDGFD_EC', run = 'martin1')

# chcek if bait in data
det <- detect(infile, 'PDGFD'); any(det) #

# Write
if (do.write){
  data <- prepare(bait = c('EC', 'PDGF'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip('Whitehead', martin = '1', bait = bait, bait.vs = 'Mock', cell=cell, note = 'V5')), quote = F, sep = '\t')
}


###########################################################
# date: 13-Jan-2019
# author: Frederik Heymann Lassen
# bait: CXCL12-beta
# control: PDGFD
# cell: THP
###########################################################
infile = "/Users/flassen/Projects/03_MICOM/data/raw/martin1//Martin_THP-PDGFD_proteins.csv"
head(read.csv(infile))
bait = 'PDGFD'
cell = 'THP'

# create directories for storaring data
dirs <- mkdir(name = 'PDGFD_THP', run = 'martin1')

# chcek if bait in data
det <- detect(infile, 'PDGFD'); any(det) #

# Write
if (do.write){
  data <- prepare(bait = c('THP', 'PDGF'), infile = infile, verbose = T, raw = T)$data
  #data <- data[,c('Accession','uniprot','rep1','rep2','imputed')]
  #colnames(data) = c('gene','accession','rep1','rep2','imputed')
  write.table(data, paste0('data/genoppi_input/',write.ip('Whitehead', martin = '1', bait = bait, bait.vs = 'Mock', cell=cell, note = 'V5')), quote = F, sep = '\t')
}




