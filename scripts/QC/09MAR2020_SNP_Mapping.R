
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
  library(igraph)
  
  #library(rProteomics, lib.loc = '~/Toolbox/rlib/')
  devtools::load_all('~/Toolbox/packages/pRoteomics/')
}

# Get loci
cad_gwas = read.table('~/Projects/03_MICOM/genelist/Nelsen2017_CAD_GWAS_Supplementary_table_4.txt', sep = '\t', header = T)
cad_gwas_soft = cad_gwas[!is.na(cad_gwas$SOFT.Pvalue),c(1,2,3,5,6,7)]
# seperate chromosome, position and whether data comes from exome seq
cad_gwas_soft$exomechip = grepl('\\*',cad_gwas_soft$Markername)
cad_gwas_soft$chr = unlist(lapply(strsplit(as.character(cad_gwas_soft$Chr.Pos), '\\:'), function(x)x[1]))
cad_gwas_soft$pos = unlist(lapply(strsplit(as.character(cad_gwas_soft$Chr.Pos), '\\:'), function(x)x[2]))
cad_gwas_soft$Chr.Pos <- NULL

# SNP to gene mapping
mapping = snp_to_gene(cad_gwas_soft[!cad_gwas_soft$exomechip,]$Markername)
mapping = data.frame(gene=names(mapping), Markername=as.character(unlist(mapping)))
mg = merge(cad_gwas_soft, mapping, by = 'Markername')


write.table(mg, 'genelist/Nelsen2017_GWAS_snp_to_genes.tsv', sep = '\t', quote = F)






