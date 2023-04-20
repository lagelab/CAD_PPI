setwd('~/Projects/03_MICOM/')
library(dplyr)
devtools::load_all('~/Projects/08_genoppi_lagelab/Genoppi/')

# setup paths
targetdate = '20NOV20'
targetrun = 'run056'
#prefix = paste0(targetdate,'_',writerun,'_micom')
paths = list.files(paste0('data/genoppi_input/', targetrun), full.names = T)
paths.whitehead = !grepl('Broad', paths)
paths.broad = grepl('Broad', paths)

# get QC'ed refs 
ref <- as.vector(read.table('run056_filtered_tier1_paths.tsv')$x)
#ref <- gsub('FLT1','FLT', ref) # they must match
#ref <- gsub('BCAS3','BCAS', ref) # they must match
#ok <- unlist(lapply(ref, function(x) x %in% basename(paths)))
ok <- unlist(lapply(basename(paths), function(x) x %in% basename(ref)))

# any missing?
if (sum(ok) != 20) {
  (bool_missing <- unlist(lapply(ref, function(x) !any(grepl(x, paths)))))
  (missing <- ref[bool_missing])  
  print(missing)
  stop('missing above IPs..') 
}

# helper funcs
check_bait <- function(df, bait){
  if (bait %in% df$gene[df$significant]){return('Yes')}
  if (bait %in% df$gene){return('No')}
  if (bait %nin% df$gene[df$significant]){return('Not found') }
}


# stats
#ndfs <- lapply(ndfs, function(x) x %>% calc_mod_ttest() %>% id_enriched_proteins())
#(median_FDR <- median(unlist(lapply(ndfs, function(x) x$FDR[x$FDR <= 0.1]))))
#(median_r <- median(unlist(lapply(ndfs, function(x) cor(x$rep1, x$rep2)))))

# genelist
haarst2017 <- data.frame(gene = as.vector(read.table('genelist/27AUG20_haarst2017_protein_coding_genelist.tsv', sep = '\t', header = T)), significant = T)
roselli <- data.frame(gene = as.vector(read.csv('data/micom_roselli.txt', header = F)), significant = T)
gnomad <- gnomad_table
gnomad$significant <- gnomad$pLI >= 0.9
gnomad$significant[is.na(gnomad$significant)] <- FALSE

summary_stats <- list()

## generate cell by cell
for (path in paths[ok]){
  
  # get summary stats 
  cat(paste0('\nReading ', path, '.. (' ,Sys.time() ,')\n'))
  bait = unlist(lapply(basename(path), function(x) unlist(strsplit(x,'\\.|vs'))[3]))
  cell = unlist(lapply(basename(path), function(x) unlist(strsplit(x,'\\.|vs'))[5]))
  facility = ifelse(grepl('Broad',unlist(lapply(basename(path), function(x) unlist(strsplit(x,'\\.|vs'))[1]))), 'Broad','Whitehead')
  dfstart = read.csv(path, sep = '\t') 
  summary_stats[[path]] <- list()
  
  # get synonyms for inweb 
  #bait_inweb = gsub('FLT', 'FLT1', bait)
  #bait_inweb = gsub('ARHGEF26\\(WT\\)', 'ARHGEF26', bait_inweb)
  #bait_inweb = gsub('BCAS', 'BCAS3', bait_inweb)
  bait_inweb = gsub('JCAD', 'KIAA1462', bait)
  
  # replace synonyms in data
  dfstart$gene[dfstart$gene == 'JCAD'] <- 'KIAA1462'
  
  for (sig.index in 1:3){
    
    # (1) FDR < 0.1
    if (sig.index == 1) {
      df = dfstart %>% id_enriched_proteins(logfc_dir = 'both', fdr_cutoff = 0.1)
      sig.name = 'LogFC.both'
    }
    
    # (2) FDR < 0.1, LogFC > 0
    if (sig.index == 2) {
      df = dfstart %>% id_enriched_proteins(logfc_dir = 'positive', fdr_cutoff = 0.1)
      sig.name = 'LogFC.positive'
    }
    
    # (3) FDR < 0.1, LogFC < 0
    if (sig.index == 3) {
      df = dfstart %>% id_enriched_proteins(logfc_dir = 'negative', fdr_cutoff = 0.1)
      sig.name = 'LogFC.negative'
    }
    
    # prepare out-file paths
    outfile = tools::file_path_sans_ext(path)
    outstatfile = paste0(outfile,'.', sig.name ,'.Enrichment.xlsx')
    outplotfile = paste0(outfile,'.', sig.name,'.Figures.pdf')
    
    # get genelist for overlap calculation
    list_inweb = get_inweb_list(bait_inweb)
    list_gnomad = gnomad
    list_haarst = haarst2017
    list_roselli = roselli

    # calculate statistics
    inweb_stats = calc_hyper(df, list_inweb, data.frame(intersectN = T), bait = bait)
    gnomad_stats = calc_hyper(df, list_gnomad, data.frame(intersectN = T), bait = bait)
    #haarst_stats = calc_adjusted_enrichment(df, list_haarst, intersectN = T, bait = bait)
    #roselli_stats = calc_adjusted_enrichment(df, list_roselli, intersectN = T, bait = bait)
    haarst_stats = calc_hyper(df, list_haarst, data.frame(intersectN = F), bait = bait)
    roselli_stats = calc_hyper(df, list_roselli, data.frame(intersectN = F), bait = bait)
    
    # combine data
    inweb_stats$statistics$successInSampleGenes = paste(inweb_stats$genes$mylist$successInSample_genes, collapse = '; ')
    gnomad_stats$statistics$successInSampleGenes = paste(gnomad_stats$genes$mylist$successInSample_genes, collapse = '; ')
    haarst_stats$statistics$successInSampleGenes = paste(haarst_stats$genes$mylist$successInSample_genes, collapse = '; ')
    roselli_stats$statistics$successInSampleGenes = paste(roselli_stats$genes$mylist$successInSample_genes, collapse = '; ')
    
    inweb_stats$statistics$list_name <- 'InWeb'
    gnomad_stats$statistics$list_name <- 'gnomAD'
    haarst_stats$statistics$list_name <- 'haarst2017'
    roselli_stats$statistics$list_name <- 'roselli'  
    
    genelist = rbind(inweb_stats$statistics, gnomad_stats$statistics, haarst_stats$statistics, roselli_stats$statistics)
    colnames(genelist) <- c("Enrichment list", "successInSample_count", "sample_count", "notSample_count", "success_count", "population_count", "p-value", "successInSampleGenes" )
    
    
    # calculate correlation
    cor12 = round(cor(df$rep1, df$rep2), 3)
    
    # write to summary stats
    summary_stats[[path]][[sig.name]] <- data.frame(
      path = basename(path),
      analysis.direction = sig.name,
      bait = bait_inweb,
      cell = cell,
      facility = facility,
      
      rep.correlation = paste(cor12),
      rep.correlation.mean = paste(round(mean(c(cor12)),3)),
      
      proteins.detected = nrow(df),
      proteins.significant = sum(df$significant),
      bait.enriched = check_bait(df, bait_inweb),
      
      inweb.detected = inweb_stats$statistics$sample_count,
      inweb.significant = inweb_stats$statistics$successInSample_count,
      inweb.pvalue = inweb_stats$statistics$pvalue,
      
      gnomad.detected = gnomad_stats$statistics$sample_count,
      gnomad.significant = gnomad_stats$statistics$successInSample_count,
      gnomad.pvalue = gnomad_stats$statistics$pvalue,
      
      haarst.detected = haarst_stats$statistics$sample_count,
      haarst.significant = haarst_stats$statistics$successInSample_count,    
      haarst.pvalue = haarst_stats$statistics$pvalue,
      
      roselli.detected = roselli_stats$statistics$sample_count,
      roselli.significant = roselli_stats$statistics$successInSample_count,   
      roselli.pvalue = roselli_stats$statistics$pvalue
      
    )
  
  }
  
}


mydf = as.data.frame(do.call(rbind, lapply(summary_stats, function(x) do.call(cbind, x))))
mydf = mydf[with(mydf, order(LogFC.positive.cell, LogFC.negative.facility)), ]
write.csv(mydf, file = '20NOV20_MICOM_Cell1_vs_Cell2_stats.csv')
write.csv(mydf, file = '~/Desktop/20NOV20_MICOM_Cell1_vs_Cell2_stats.csv')


# VENN Diagrams
####################
# Helper functions #
####################


# collapse a list of IPs into a superset
collapse_ips <- function(lst){
  superset = as.data.frame(do.call(rbind, lapply(lst, function(df){df[,c('gene','significant')]})))
  sig_genes = unique(superset$gene[superset$significant])
  superset$both_insig_and_insig = superset$gene %in% sig_genes & superset$significant == FALSE
  superset$significant = superset$significant | superset$both_insig_and_insig
  outdf = superset[,c('gene','significant')]
  outdf = outdf[!duplicated(outdf), ]
  return(outdf)
}


# setup baits, cells, and data

# data
data_igg <- lapply(paths[ok], function(x) read.csv(x, sep = '\t'))
thenames <- basename(paths)[ok]

# get all baits and cells
all_baits <- unlist(lapply(thenames, function(x) unlist(strsplit(x, '\\.|vs'))[3]))
all_baits <- gsub('^BCAS$','BCAS3', all_baits)
all_baits <- gsub('^FLT$','FLT1', all_baits) # they must match
all_baits <- gsub('ARHGEF26\\(WT\\)','ARHGEF26', all_baits) 
all_baits <- gsub('JCAD','KIAA1462', all_baits) 
all_cells <- unlist(lapply(thenames, function(x) unlist(strsplit(x, '\\.'))[4]))
newnames <- paste(all_baits, all_cells, sep = '.')
names(data_igg) <- newnames

# dont repeat any when we loop
all_cells <- unique(all_cells)
all_baits <- unique(all_baits)

mysummary <- list()

for (bait in all_baits){
  
  inweb_bait <- genoppi::get_inweb_list(bait)
  
  #####################
  ### cell versus cell
  
  for (cell1 in all_cells){
    for (cell2 in all_cells){
      
      # setup regex
      RE_bait = grepl(bait, names(data_igg))
      RE_cell1 = grepl(cell1, names(data_igg))
      RE_cell2 = grepl(cell2, names(data_igg))
      
      # get paths
      path_bait_cell1 = names(data_igg)[RE_bait & RE_cell1]
      path_bait_cell2 = names(data_igg)[RE_bait & RE_cell2]
      
      # only continue, if the bait is represented in each cell
      if (!identical(path_bait_cell1, character(0)) & !identical(path_bait_cell2, character(0))){
        
        # get data
        data_bait_cell1 = lapply(path_bait_cell1, function(p) data_igg[[p]] %>% id_enriched_proteins())
        data_bait_cell2 = lapply(path_bait_cell2, function(p) data_igg[[p]] %>% id_enriched_proteins())
        
        # collapse to superset
        superset_bait_cell1 = collapse_ips(data_bait_cell1)
        superset_bait_cell2 = collapse_ips(data_bait_cell2)
        
        # calculate statistics
        index = paste0(bait,':',cell1,'vs',cell2)
        hypergeom = calc_hyper(superset_bait_cell1, superset_bait_cell2, data.frame(listName = index,intersectN = T), bait = bait)
        
        # GTEX
        
        
        mysummary[[index]] <- hypergeom$statistics
        
        
        
      }
      
    }
  }
  
  
  ## cell versus inweb
  for (cell in all_cells){
    
    # setup regex
    RE_bait = grepl(bait, names(data_igg))
    RE_cell = grepl(cell, names(data_igg))
    
    # get paths
    path_bait_cell = names(data_igg)[RE_bait & RE_cell]
    
    if (!identical(path_bait_cell, character(0))){
      
      # get data
      data_bait_cell = lapply(path_bait_cell, function(p) data_igg[[p]] %>% id_enriched_proteins())
      
      # collapse to superset
      superset_bait_cell = collapse_ips(data_bait_cell)
      
      # calculate statistics
      index = paste0(bait,':',cell,'vsInWeb_InBiomap')
      hypergeom = calc_hyper(superset_bait_cell, inweb_bait, data.frame(listName = index,intersectN = T), bait = bait)
      mysummary[[index]] <- hypergeom$statistics
      
    }
    
  }
  
  
  
  ##########################
  ### InWeb vs all cell types
  
  
  # setup regex
  RE_bait = grepl(bait, names(data_igg))
  
  # get paths
  path_bait_cell = names(data_igg)[RE_bait]
  
  # get data
  data_bait_cell = lapply(path_bait_cell, function(p) data_igg[[p]] %>% id_enriched_proteins())
  
  # collapse to superset
  superset_bait_cell = collapse_ips(data_bait_cell)
  
  # calculate statistics
  index=paste0(bait,':',paste0(all_cells,collapse='+'), 'vsInWeb_InBiomap')
  hypergeom = calc_hyper(superset_bait_cell, inweb_bait, data.frame(listName = index, intersectN = T), bait = bait)
  mysummary[[index]] <- hypergeom$statistics
  
}


df = data.frame(do.call(rbind, mysummary))
df = data.frame(apply(df, 2, function(x) unlist(x)))
#df$celltype <- ifelse(df$area1 == 'HDFN' & df$area2 == 'InWeb_InBiomap', 'neuron', NA)
#df$celltype <- ifelse(grepl("(A375)|(G401)|(T47D)", df$list_name), 'cancer', 'neuron')
df$bait <- unlist(lapply(strsplit(df$list_name, '\\:'), function(x) x[1]))
df$group1 <- unlist(lapply(strsplit(df$list_name, '(\\:)|(vs)'), function(x) x[2]))
df$group2 <- unlist(lapply(strsplit(df$list_name, '(\\:)|(vs)'), function(x) x[3]))
df$pvalue.lt.0.05 <- ifelse(as.numeric(df$pvalue) < 0.05, 'Y','N')
write.csv(df, file = paste0(targetdate,'_micom_overlap_table.csv'), quote = F)



## helpers
compare <- function(what, xname, yname, breaks = 20){
  x = data_igg[[xname]]
  y = data_igg[[yname]]
  #x = x[x$significant, ]
  #y = y[y$significant, ]
  lims = range(unlist(list(x[[what]], y[[what]])))*1.2
  p1 <- hist(x[[what]], breaks = breaks, xlim = lims,  main = xname, xlab = what, col = rgb(0, 0, 1, 1/4))
  p2 <- hist(y[[what]], breaks = breaks, xlim = lims, main = yname, xlab = what, col = rgb(1, 0, 0, 1/4))
  lims_count <- range(c(p1$counts, p2$counts))
  plot(p1, col = rgb(0, 0, 1, 1/4), xlim = lims, xlab = what, main = paste(xname, 'vs', yname))
  plot(p2, col = rgb(1, 0, 0, 1/4), xlim = lims, add = T)
}

compare.logFC <- function(xname, yname){
  x = data_igg[[xname]][,c('gene', 'logFC')]
  y = data_igg[[yname]][,c('gene', 'logFC')]
  res = merge(x, y, by = 'gene')
  r = cor(res$logFC.x, res$logFC.y)
  plot(x = res$logFC.x, y = res$logFC.y, xlab = paste(xname,'logFC'), ylab = paste(yname,'logFC'), main = paste(xname,'vs',yname ,'r =', round(r,4)))
}


## Compare LogFC beteween IPs when hypergeometric overlap for single bait in differnt cell type is significant.
names(data_igg)

thebaits = c('FN1', 'PHACTR1', 'PLPP3', 'ARHGEF26', 'HDAC9', 'EDN1')

pdf('20NOV20_logfc_comparison_of_selcted_baits_only_significant.pdf', width = 9, height = 14)
par(mfrow=c(6,3))
for (b in thebaits){
  
  xname = paste0(b, '.EC')
  yname = paste0(b, '.SMC')
  compare('logFC', xname, yname)
  
}
graphics.off()

pdf('20NOV20_logfc_rep_comparison.pdf', width = 10, height = 7)

par(mfrow=c(2,3))
for (b in thebaits){
  
  xname = paste0(b, '.EC')
  yname = paste0(b, '.SMC')
  compare.logFC(xname, yname)
  
}
graphics.off()










