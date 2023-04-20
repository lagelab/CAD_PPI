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
list.files('data/raw')


######################
# KCNK (EC) ANALYSIS #
######################
# Martin 3

# ready data and baits
infile = 'data/raw/martin3/3martin_BCAS3_proteins.csv'
bait = 'BCAS3'
cell = 'EC'
impute = 'SHIFT'
  
# prepare data and ready directories
dirs <- mkdir(bait, cell) 
  
# Prepare data for analysis by getting replicate folds
df <- prepare(c(cell, bait), infile = infile, impute = list(shift = -1.8, stdwidth = 0.3), verbose = T)
colnames(df)[1] <- 'gene'
  
# get logfoldchange and do moderated ttest
data <- df %>% mttest()
write.table(data,file=dirs$txtpath,row.names=F,sep="\t",quote=F)
#pdf(dirs$pdfpath, height=4, width=4)
  
## plot volcano and scatter plot
data %>% designate(FDR < 0.1) %>% plotVolcano(bait, title = 'BCAS3 vs Control (ES)')
data %>% designate(FDR < 0.1) %>% plotScatter(bait, title = 'BCAS3 Replicate Correlation')

# plot imputed data
data %>% designate(imputed == TRUE, FDR < 0.1) %>% plotVolcano(bait, title = 'Imputation summary', sub2 ='imputed')
  
### plot inweb interactors
known.interactors = interactors('BCAS3', T)
data %>% designate(FDR < 0.1) %>% plotOverlap(bait, known.interactors, 'BCAS3 vs Control (ES): InWeb Overlap')
graphics.off()
  
######################
# KCNK5 (EC) ANALYSIS #
######################
# Martin 3

# ready data and baits
infile = 'data/raw/martin3/3martin_KCNK5_proteins.csv'
bait = 'KCNK5'
cell = 'EC'
impute = 'SHIFT'

# prepare data and ready directories
dirs <- mkdir(bait, cell, impute) 

# Prepare data for analysis by getting replicate folds
df <- prepare(c(cell, bait), infile = infile, impute = list(shift = -1.8, stdwidth = 0.3), verbose = T)
colnames(df)[1] <- 'gene'

# get logfoldchange and do moderated ttest
data <- df %>% mttest()
write.table(data,file=dirs$txtpath,row.names=F,sep="\t",quote=F)
pdf(dirs$pdfpath, height=4, width=4)

## plot volcano and scatter plot
data %>% designate(FDR < 0.1) %>% plotVolcano(bait, title = 'KCNK5 vs Control (ES)')
data %>% designate(FDR < 0.1) %>% plotScatter(bait, title = 'KCNK5 Replicate Correlation')

# plot imputed data
data %>% designate(imputed == TRUE, FDR < 0.1) %>% plotVolcano(bait, title = 'Imputation summary', sub2 ='imputed')

### plot inweb interactors
known.interactors = interactors('KCNK5', T)
data %>% designate(FDR < 0.1) %>% plotOverlap(bait, known.interactors, 'KCNK5 vs Control (ES): InWeb Overlap')
graphics.off()

######################
# KSR2 (EC) ANALYSIS #
######################
# martin 4

# ready data and baits
infile = 'data/raw/martin4/4martin_KSR2_proteins.csv'
bait = 'KSR'
cell = 'EC'
impute = 'SHIFT'

# prepare data and ready directories
dirs <- mkdir(name='EC_KSR2') 

# Prepare data for analysis by getting replicate folds
df <- prepare(c(cell, bait), infile = infile, impute = list(shift = -1.8, stdwidth = 0.3), verbose = T)
colnames(df)[1] <- 'gene'

# get logfoldchange and do moderated ttest
data <- df %>% mttest()
write.table(data,file=dirs$txtpath,row.names=F,sep="\t",quote=F)
pdf(dirs$pdfpath, height=4, width=4)

## plot volcano and scatter plot
data %>% designate(FDR < 0.1) %>% plotVolcano(bait, title = 'KSR2 vs Control (ES)')
data %>% designate(FDR < 0.1) %>% plotScatter(bait, title = 'KSR2 replicate correlations')

# plot imputed data
data %>% designate(imputed == TRUE, FDR < 0.1) %>% plotVolcano(bait, title = 'Imputation summary', sub2 ='imputed')

### plot inweb interactors
known.interactors = interactors('KSR2', T)
data %>% designate(FDR < 0.3) %>% plotOverlap(bait, known.interactors, 'KSR2 vs Control (ES): InWeb Overlap')
graphics.off()

#######################
# ARHGEF26(WT-mutant) #
#######################
# martin 4

### wild type

# ready data and baits
infile = 'data/raw/martin4/4martin_ARHGEF26(WT-mutant)_proteins.csv'
bait = 'WT'
cell = 'EC'
impute = 'SHIFT'

# prepare data and ready directories
dirs <- mkdir(name = 'EC_ARHGEF26(WT)') 

# Prepare data for analysis by getting replicate folds
df <- prepare(c(cell, bait), infile = infile, impute = list(shift = -1.8, stdwidth = 0.3), verbose = T)
colnames(df)[1] <- 'gene'

# get logfoldchange and do moderated ttest
data <- df %>% mttest()
write.table(data,file=dirs$txtpath,row.names=F,sep="\t",quote=F)
pdf(dirs$pdfpath, height=4, width=4)

## plot volcano and scatter plot
data %>% designate(FDR < 0.1) %>% plotVolcano(bait, title = 'ARHGEF26(WT) vs Control (ES)')
data %>% designate(FDR < 0.1) %>% plotScatter(bait, paste(bait, 'ARHGEF26(WT) Replicate Correlations'))

# plot imputed data
data %>% designate(imputed == TRUE, FDR < 0.1) %>% plotVolcano(bait, title = 'Imputation summary', sub2 ='imputed')

### plot inweb interactors
known.interactors = interactors('ARHGEF26', T)
data %>% designate(FDR < 0.1) %>% plotOverlap(bait, known.interactors, 'ARHGEF26(WT) vs Control (ES): InWeb Overlap')
graphics.off()


### Mutant

# ready data and baits
infile = 'data/raw/martin4/4martin_ARHGEF26(WT-mutant)_proteins.csv'
bait = 'mt'
cell = 'EC'
impute = 'SHIFT'

# prepare data and ready directories
dirs <- mkdir(name = 'EC_ARHGEF26(MT)') 

# Prepare data for analysis by getting replicate folds
df <- prepare(c(cell, bait), infile = infile, impute = list(shift = -1.8, stdwidth = 0.3), verbose = T)
colnames(df)[1] <- 'gene'

# get logfoldchange and do moderated ttest
data <- df %>% mttest()
write.table(data,file=dirs$txtpath,row.names=F,sep="\t",quote=F)
pdf(dirs$pdfpath, height=4, width=4)

## plot volcano and scatter plot
data %>% designate(FDR < 0.1) %>% plotVolcano(bait, title = 'ARHGEF26(MT) vs Control (ES)')
data %>% designate(FDR < 0.1) %>% plotScatter(bait, paste(bait, 'ARHGEF26(MT) Replicate Correlations'))

# plot imputed data
data %>% designate(imputed == TRUE, FDR < 0.1) %>% plotVolcano(bait, title = 'Imputation summary', sub2 ='imputed')

### plot inweb interactors
known.interactors = interactors('ARHGEF26', T)
data %>% designate(FDR < 0.1) %>% plotOverlap(bait, known.interactors, 'ARHGEF26(WT) vs Control (ES): InWeb Overlap')
graphics.off()





#'data/raw/martin4/4martin_ARHGEF26(WT-mutant)_proteins.csv',
#'data/raw/martin4/4martin_ARHGEF26_proteins.csv',





