#.libPaths('~/Toolbox/rlib')
setwd('~/Projects/03_MICOM/')
devtools::load_all('~/Projects/08_genoppi_lagelab/Genoppi/')
paths = list.files(paste0('data/genoppi_input/', 'run019'), full.names = T)
paths[grepl('HEK',paths)]



## coding variant differential expression in HEK cells
fpath <-


#fpath = "data/genoppi_input/run019/Genoppi_Whitehead_martin9_MinImputed.ARHGEF26(WT)vsARHGEF26(MT).SMC.3xFlag.21JUL2020.tsv" 
fpath = "data/genoppi_input/run019/Genoppi_Whitehead_martin9_MinImputed.ARHGEF26(WT)vsARHGEF26(MT).teloHAEC.3xFlag.21JUL2020.tsv"
df <- read.csv(fpath, sep = '\t') %>% id_enriched_proteins(logfc_dir = 'both')
df[df$gene == 'ARHGEF26',]


plot_volcano_basic(df) %>%
  plot_overlay(as.bait('ARHGEF26'))

