## initial analysis
.libPaths('~/Toolbox/rlib')
setwd('~/Projects/03_MICOM/')
library(dplyr)
library(limma)
library(ggplot2)
library(ggrepel)

#library(rProteomics, lib.loc = '~/Toolbox/rlib/')
devtools::load_all('~/Toolbox/packages/pRoteomics/')

## run analysis
vbait <- c('BCAS3',
          'KCNK5')
vcell <- c('EC',
          'EC')
vimpute <- c('SHIFT',
            'SHIFT')
#dat_bcas3_pp <- 'data/raw/martin3/3martin_BCAS3_protein-peptides.csv'
vdat <- c('data/raw/martin3/3martin_BCAS3_proteins.csv',
         'data/raw/martin3/3martin_KCNK5_proteins.csv')



for (i in 1:2){
  
  # prepare data and get replicates
  bait = vbait[i]
  cell = vcell[i]
  impute = vimpute[i]
  dat = vdat[i]
  
  df <- prepare(c(cell, bait), infile = dat, impute = list(shift = -0.8, stdwidth = 0.3))
  
  # generate directories
  newdir <- paste0(c('derived/',bait,'_', cell), collapse = '')
  output <- paste0(c(bait,'_', cell,'_IMPUTE=',impute,'_',toupper(format(Sys.time(), "%d%b%Y")),'.txt'), collapse = '')
  newfile <- file.path(newdir, output)
  dir.create(newdir)
  colnames(df) = c('gene', 'imputed', 'rep1', 'rep2')
  write.table(df, newfile, sep = '\t')
  
  
  # genoppi standardized pipeline
  output <- paste0(c(bait,'_', cell,'_IMPUTE=',impute,'_',toupper(format(Sys.time(), "%d%b%Y"))), collapse = '')
  genoppi(newfile, bait_name = bait, debug = F, p_cutoff = NA, fdr_cutoff = 0.1,
          output_stats_file = file.path(newdir, paste0(output,'_statistics.txt')),
          output_plots_file = file.path(newdir, paste0(output,'_plots.pdf')))
  
}







