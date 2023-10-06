# get all relevant data files

#m15 <- list.files('~/Projects/03_MICOM/data/genoppi_input/raw/martin9/', full.names = T)
#files <- c(m9, m15)
#files <- files[grepl('mut', files) & grepl('proteins\\.csv', files)]
paths <- list.files(paste0('data/genoppi_input/', 'run019'), full.names = T)
files <- paths[grepl('ARHGEF26\\((W|M)T\\)vsMock',paths)]

xp = lapply(files, function(f) {
  x = read.csv(f, sep = '\t')
  return(cor(x$rep1, x$rep2))
})

names(xp) <- basename(files)

## SMC 

f1 = "data/genoppi_input/run019/Genoppi_Whitehead_martin9_MinImputed.ARHGEF26(MT)vsMock.SMC.3xFlag.21JUL2020.tsv" 
f2 = "data/genoppi_input/run019/Genoppi_Whitehead_martin9_MinImputed.ARHGEF26(WT)vsMock.SMC.3xFlag.21JUL2020.tsv" 
f3 = "data/genoppi_input/run019/Genoppi_Whitehead_martin9_MinImputed.ARHGEF26(WT)vsARHGEF26(MT).SMC.3xFlag.21JUL2020.tsv" 

mt = read.csv(f1, sep = '\t') %>% id_enriched_proteins(logfc_dir = 'positive')
wt = read.csv(f2, sep = '\t') %>% id_enriched_proteins(logfc_dir = 'positive')
mtwt = read.csv(f3, sep = '\t') %>% id_enriched_proteins(logfc_dir = 'both')

prot_mt <- list(highlight = data.frame(gene = mt$gene[mt$significant & mt$gene %nin% 'ARHGEF26'], significant = TRUE, col_significant = "#41AB5D"))
prot_wt <- list(highlight = data.frame(gene = wt$gene[wt$significant & wt$gene %nin% 'ARHGEF26'], significant = TRUE, col_significant = "#41AB5D"))
prot_mtwt <- list(highlight = data.frame(gene = mtwt$gene[mtwt$significant], significant = TRUE, col_significant = "#41AB5D"))

mt %>% plot_scatter_basic() %>% plot_overlay(as.bait('ARHGEF26')) %>% plot_overlay(prot_mt)  #%>% make_interactive()
wt %>% plot_scatter_basic() %>% plot_overlay(as.bait('ARHGEF26')) %>% plot_overlay(prot_wt) #%>% make_interactive()
mtwt %>% plot_scatter_basic() %>% plot_overlay(as.bait('ARHGEF26')) %>% plot_overlay(prot_mtwt) # %>% make_interactive()


goi = 'ALYREF' # observe gene of intersest #MDK
mt %>% plot_volcano_basic() %>% plot_overlay(as.bait('ARHGEF26')) %>% plot_overlay(as.goi(goi)) %>% make_interactive()
wt %>% plot_volcano_basic() %>% plot_overlay(as.bait('ARHGEF26')) %>% plot_overlay(as.goi(goi)) %>% make_interactive()
mtwt %>% plot_volcano_basic() %>% plot_overlay(as.bait('ARHGEF26'))  %>% make_interactive()

# what proteins are significant in mutant
mt$gene[mt$significant]
wt$gene[wt$significant]
mtwt$gene[mtwt$significant == T & mtwt$logFC > 0]

# explore original data file
list.files('~/Projects/03_MICOM/data/genoppi_input/raw/martin15/', full.names = T)
p <- "/Users/flassen/Projects/03_MICOM/data/genoppi_input/raw/martin15//15martin_ARHGEF(WT&mut)_NoNorm_First_proteins.csv"
#p <- "/Users/flassen/Projects/03_MICOM/data/genoppi_input/raw/martin9//9martin_ARHRGF(WTmut)-HAEC_protein-peptides.csv"  
#p <- "/Users/flassen/Projects/03_MICOM/data/genoppi_input/raw/martin9//9martin_ARHRGF(STmut)-teloHAEC_proteins.csv" 
f <-read.csv(p)
f[grepl('Martin', f$Description),]
f[grepl('FLAG', f$Accession),]
f[grepl('ARHGEF26', f$Accession),]

## explore HEK cells
myfiles <- list.files('~/Projects/03_MICOM/data/genoppi_input/raw/', full.names = T, recursive = T)
myfiles[grepl('HEK',myfiles)]

p <- "/Users/flassen/Projects/03_MICOM/data/genoppi_input/raw//martinBroad/ARHGEF26/ARHGEF26_HEK293T_nonSGS.csv" 
f <-read.csv(p, sep = '\t')
f[grepl('ARHGEF',f$gene),]

