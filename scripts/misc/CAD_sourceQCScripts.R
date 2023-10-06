## initial analysis
.libPaths('~/Toolbox/rlib')
setwd('~/Projects/03_MICOM/')
library(dplyr)
library(limma)
library(ggplot2)
library(ggrepel)
library(hash)
library(dplyr)
library(openxlsx)
  
#library(rProteomics, lib.loc = '~/Toolbox/rlib/')
devtools::load_all('~/Toolbox/packages/pRoteomics/')

## load CAD geneset
geneset <- read.table('~/Toolbox/datasets/CAD/cad_gwas_genes_prio1_2.txt');
geneset <- data.frame(genes=geneset, significant=TRUE)
colnames(geneset) <- c('gene','significant')

do.write = T

warn = function(x) write(x, stderr())

write.ip = function(martin, bait, cell, bait.vs = 'Mock', facility = 'Whitehead', note = 'expanded',imputation = 'MinImputed'){
  return(paste0('Genoppi_',facility,'_martin',martin,'_',imputation,'.',bait,'vs',bait.vs,'.',cell,'.',note,'.tsv'))
}

files <- list.files('~/Projects/03_MICOM/scripts/process_intensity_data/', pattern = 'martin', full.names = T)
for (f in files){
  cat('\n')
  cat(f)
  cat('\n')
  
  source(f)
}


# check result and rename
if (F){
  
  mypaths = '~/Projects/03_MICOM/data/genoppi_input/run050'
  newfiles <- list.files(mypaths)
  #unlist(lapply(strsplit(newfiles, '\\.'), function(x) x[2]))
  
  # rename
  newnames <- gsub('Genoppi_Whitehead_mar', 'm',newfiles)
  newnames <- gsub('\\_MinImputed', '',newnames)
  newnames <- gsub('mtinBroad', 'Broad\\.mB',newnames)
  newnames <- gsub('mtin', 'Whitehead\\.m',newnames)
  
  ndf = data.frame(old = newfiles, new = newnames)

  for (i in 1:nrow(ndf)){
    
    ffrom = paste0(mypaths, '/', ndf$old[i])
    fto = paste0(mypaths, '/', ndf$new[i])
    df = read.csv(ffrom, sep = '\t')
    write.table(df, fto, sep = '\t')
    write(paste0(fto,'... [W]'), stdout())
    unlink(ffrom)
    write(paste0(ffrom,'... [D]'), stdout())
    
  }


}