setwd('~/Projects/03_MICOM/')

devtools::load_all('~/Projects/08_genoppi_lagelab/Genoppi/')

# what proteins do we have here
bp <- goa_bp_table[,c(1,3)]
colnames(bp) <- c('gene', 'dataset')
bp$significant <- FALSE

mf <- goa_mf_table[,c(1,3)]
colnames(mf) <- c('gene', 'dataset')
mf$significant <- FALSE

c2 <- msigdb_c2_table
colnames(c2) <- c('gene', 'dataset')
c2$significant <- FALSE

analysis <- list(bp = bp, mf = mf, c2 = c2)

# setup paths and run data
paths = read.csv('run056_filtered_tier1_paths.tsv', sep = '\t')$x
target = 'enrichment_analysis/run056/'

for (i in seq_along(analysis)){

  outdir = file.path(target,names(analysis)[i])
  dir.create(outdir)
  
  for (p in paths){
    
    # prepare names for outfile
    name = basename(p)
    df = read.csv(p, sep = '\t')
    bait = unlist(strsplit(name, '\\.|vs'))[3]
    cell = unlist(strsplit(name, '\\.|vs'))[5]
    newname = gsub('20NOV2020\\.tsv', paste0('count_enrichment_',names(analysis)[i], '.tsv'), name)
    outfile = file.path(outdir,newname)
    
    if (!file.exists(outfile)){
      res = apply_calc_hyper(df, analysis[[i]], col.by = 'dataset', bait = bait, verbose = T)
      write.table(res, outfile, sep = '\t', row.names = F, quote = F)
    }

  }

}





