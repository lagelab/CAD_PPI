
if (F){
  ## initial analysis
  .libPaths('~/Toolbox/rlib')
  setwd('~/Projects/03_MICOM/')
  devtools::load_all('~/Toolbox/packages/pRoteomics/')
  devtools::load_all('~/Projects/08_genoppi_lagelab/Genoppi/')
  library(dplyr)
}


######## data directory
(files <- list.files('~/Projects/03_MICOM/data/raw/martin17', pattern = 'proteins\\.csv', full.names = T))
stopifnot(length(files) == 7)





### ADAMTS7 (EC)
(infile = files[1])
file = read.csv(infile)
colnames(file)

# extract replicates
bait = 'ADAMTS7'
cell = 'EC'
dirs <- mkdir(bait = bait, cell = cell, run = 'martin17')
data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, 
                cols = c('Accession', 'Intensity.EC.A7.1.iTRAQ4.114.', 'Intensity.EC.mockA1.iTRAQ4.116.',
                         'Intensity.EC.A7.2.iTRAQ4.115.', 'Intensity.EC.mockA2.iTRAQ4.117.'))$data

pairs(data[,grepl('Intensity.iTRAQ4',colnames(data))])
plot(data$rep1, data$rep2, main = 'ADAMTS7 (EC)')

write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '17', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')

### BCAS3-HAEC
(infile = files[2])
file = read.csv(infile)
colnames(file)

# extract replicates
bait = 'BCAS3'
cell = 'EC'
dirs <- mkdir(bait = bait, cell = cell, run = 'martin17')
data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, 
                cols = c('Accession', 
                         'Intensity.EC.BCAS.1.iTRAQ4.114.', 'Intensity.EC.mockB1.iTRAQ4.116.',
                         'Intensity.EC.BCAS.2.iTRAQ4.115.', 'Intensity.EC.mockB2.iTRAQ4.117.'))$data

pairs(data[,grepl('Intensity.iTRAQ4',colnames(data))])
plot(data$rep1, data$rep2, main = 'BCAS3 (EC)')

write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '17', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')

### EDNRA-HAEC
(infile = files[3])
file = read.csv(infile)
colnames(file)

# extract replicates
bait = 'EDNRA'
cell = 'EC'
dirs <- mkdir(bait = bait, cell = cell, run = 'martin17')
data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, 
                cols = c('Accession', 
                         'Intensity.EC.EDNRA1.iTRAQ4.114.', 'Intensity.EC.mockE1.iTRAQ4.116.',
                         'Intensity.EC.EDNRA2.iTRAQ4.115.', 'Intensity.EC.mockE2.iTRAQ4.117.'))$data

pairs(data[,grepl('Intensity.iTRAQ4',colnames(data))])
plot(data$rep1, data$rep2, main = 'EDNRA (EC)')

write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '17', bait = bait, cell=cell, bait.vs = 'Mock', note = 'endogenous')), quote = F, sep = '\t')


### FLT1-SMC
(infile = files[4])
file = read.csv(infile)
colnames(file)

# extract replicates
bait = 'FLT1'
cell = 'SMC'
dirs <- mkdir(bait = bait, cell = cell, run = 'martin17')
data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, 
                cols = c('Accession', 
                         'Intensity.SMC.FLT1.iTRAQ4.114.', 'Intensity.SMC.mockF1.iTRAQ4.116.',
                         'Intensity.SMC.FLT2.iTRAQ4.115.', 'Intensity.SMC.mockF2.iTRAQ4.117.'))$data

pairs(data[,grepl('Intensity.iTRAQ4',colnames(data))])
plot(data$rep1, data$rep2, main = 'FLT1 (SMC)')

write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '17', bait = bait, cell=cell, bait.vs = 'Mock', note = 'endogenous')), quote = F, sep = '\t')



### JCAD-EC
(infile = files[5])
file = read.csv(infile)
colnames(file)

# extract replicates
bait = 'JCAD'
cell = 'EC'
dirs <- mkdir(bait = bait, cell = cell, run = 'martin17')
data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, 
                cols = c('Accession', 
                         'Intensity.EC.JCAD1.iTRAQ4.114.', 'Intensity.EC.mockJ1.iTRAQ4.116.',
                         'Intensity.EC.JCAD2.iTRAQ4.115.', 'Intensity.EC.mockJ2.iTRAQ4.117.'))$data

pairs(data[,grepl('Intensity.iTRAQ4',colnames(data))])
plot(data$rep1, data$rep2, main = 'JCAD (EC)')

write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '17', bait = bait, cell=cell, bait.vs = 'Mock', note = 'endogenous')), quote = F, sep = '\t')

### KSR2-EC
(infile = files[6])
file = read.csv(infile)
colnames(file)

# extract replicates
bait = 'KSR2'
cell = 'EC'
dirs <- mkdir(bait = bait, cell = cell, run = 'martin17')
data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, 
                cols = c('Accession', 
                         'Intensity.EC.KSR.1.iTRAQ4.114.', 'Intensity.EC.mock.1.iTRAQ4.116.',
                         'Intensity.EC.KSR.2.iTRAQ4.115.', 'Intensity.EC.mock.2.iTRAQ4.117.'))$data

pairs(data[,grepl('Intensity.iTRAQ4',colnames(data))])
plot(data$rep1, data$rep2, main = 'KSR (EC)')

write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '17', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')



### KSR2-SMC
(infile = files[7])
file = read.csv(infile)
colnames(file)

# extract replicates
bait = 'KSR2'
cell = 'SMC'
dirs <- mkdir(bait = bait, cell = cell, run = 'martin17')
data <- prepare(c(cell, bait), infile = infile, verbose = T, raw = T, 
                cols = c('Accession', 
                         'Intensity.SMC.KSR.1.iTRAQ4.114.', 'Intensity.SMC.mockK.1.iTRAQ4.116.',
                         'Intensity.SMC.KSR.2.iTRAQ4.115.', 'Intensity.SMC.mockK.2.iTRAQ4.117.'))$data

pairs(data[,grepl('Intensity.iTRAQ4',colnames(data))])
plot(data$rep1, data$rep2, main = 'KSR2 (EC)')

write.table(data, paste0('data/genoppi_input/',write.ip(facility = 'Whitehead', martin = '17', bait = bait, cell=cell, bait.vs = 'Mock', note = '3xFLAG')), quote = F, sep = '\t')







