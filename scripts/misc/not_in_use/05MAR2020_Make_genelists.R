if (F){
  .libPaths('~/Toolbox/rlib')
  setwd('~/Projects/03_MICOM/')
  library(dplyr)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(hash)
  library(dplyr)
  library(igraph)
  library(data.table)
  
  #library(rProteomics, lib.loc = '~/Toolbox/rlib/')
  devtools::load_all('~/Toolbox/packages/pRoteomics/')
}

#
c4d_cad_gwas = fread('~/Projects/03_MICOM/genelist/C4D_CAD_DISCOVERY_METAANALYSIS_UPDATE.TXT', header = T)
cardiogram_gwas = fread('~/Projects/03_MICOM/genelist/CARDIoGRAM_GWAS_RESULTS.txt', header = T)

sum(c4d_cad_gwas$PVALUE < 5e-8)



