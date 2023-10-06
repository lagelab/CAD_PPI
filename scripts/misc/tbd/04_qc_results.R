setwd('~/Projects/03_MICOM/')
library(dplyr)
devtools::load_all('~/Projects/08_genoppi_lagelab/Genoppi/')

# ready paths and prefixes
targetdate = '201111'
targetrun = 'run033'
writerun = 'run034'

# setup paths
prefix = paste0(targetdate,'_',writerun,'_micom')
paths = list.files(paste0('data/genoppi_input/', targetrun), full.names = T)
paths.whitehead = !grepl('Broad', paths)
paths.broad = grepl('Broad', paths)

# check datamand remove rows that do not pass checks
dfs = lapply(paths, function(path) {
  df = read.csv(path, sep = '\t', stringsAsFactors = F)
  df = df[df$remove == FALSE,]
  return(df)
})

## QC stats

# how many genes are present
genes = unlist(lapply(dfs, function(x) x$gene.final))
n = length(genes)
print(paste(n, 'genes with', length(unique(genes)), 'unique genes')) # "65023 genes with 5280 unique genes"

# how many genes are not an approvced HGNC symbol?
hgnc.aliases = read.csv( 'data/hgnc_aliases.tsv', sep = '\t') 
(hgnc.table = sort(table(genes[!genes %in% hgnc.aliases$from]))) # 80 total and 34 unique.

# how many genes are NAs?
print(paste(sum(is.na(genes)), 'genes are not mapped!'))
lapply(dfs, function(x) x[is.na(x$gene.final), ])

# how many genomic locations are not mapped?
chrom = unlist(lapply(dfs, function(x) x[is.na(x$chrom),'gene.final']))
chrom = unlist(lapply(dfs, function(x) x$chrom))
print(paste(sum(is.na(chrom)), 'genes are are not mapped to a chromsomal position!')) # 78 (31 unique)

# how many reps are NA?
rep.na = unlist(lapply(dfs, function(x) sum(is.na(x$rep1) | is.na(x$rep2))))
print(paste(rep.na, 'data had NAs in rep..'))

# how many values are imputed?
intensity.imputed = unlist(lapply(dfs, function(x) sum(na.omit(x$imputed))))
print(paste(sum(intensity.imputed), 'intesity values were imputed!'))

# how many baits can't be found?
baits = c('BCAS3' ,'FLT1', 'KCNK5', 'ARHGEF26', 'JCAD', 'FN1', 'EDNRA', 'HDAC9', 'PLPP3', 'ADAMTS7','PHACTR1', 'KSR2')
bait.aliases = c('KIAA1462', 'PPAP2B', 'RPEL1')
(n.baits = sum(!baits %in% genes)) # all baits found in IPs (zero was not found)

# can bait HGNC aliases be found?
what = c(bait.aliases, hgnc.aliases[hgnc.aliases$from %in% baits,]$to)
(n.baits.alises = sum( what %in% genes)) # 0 aliases found found


#sum(unlist(lapply(dfs, function(x) any(colnames(x) %in% 'rep1'))))
#sum(unlist(lapply(dfs, function(x) any(colnames(x) %in% 'rep2'))))
#x <- unlist(lapply(dfs, function(x) any(colnames(x) %in% 'unique.peptides')))
#lapply(dfs, function(x) any(colnames(x) %in% 'X.unique'))

# we should not order the moderated t-test output
calc_mod_ttest_new <- function (df) 
{
  myfit <- limma::lmFit(subset(df, select = -c(gene)), method = "robust", maxit = 1000)
  myfit <- limma::eBayes(myfit)
  modtest <- limma::topTable(myfit, number = nrow(myfit), sort.by = "none")
  colnames(modtest)[4:5] <- c("pvalue", "FDR")
  result <- data.frame(cbind(df, modtest[, -c(2, 3, 6)]))
  #result <- result[with(result, order(-logFC, FDR)), ]
  return(result)
}



path_bait_table <- read.csv('data/pathname_bait_to_gene.csv', stringsAsFactors = F)
translate_path_bait <- function(x, table) {
  return(table$To[table$From == x])
}

roselli.gene <- as.vector(unlist(read.csv('data/micom_roselli.txt', header = F)))



# function to standarize bait names to inweb
standardize <- function(x){
  x = gsub('PLPP3','PPAP2B', x)
  return(x)
}




## QC every individual bait and write file
#dfs = lapply(paths, function(path) {

summary <- list()

for (path in paths){
  
  # check data
  df = read.csv(path, sep = '\t', stringsAsFactors = F) 
  path = gsub(targetrun, writerun, path)
  newpath = path
  rows.qc.total <- nrow(df)
  rows.qc.discarded = sum(df$remove == FALSE)
  
  # remove rows here..
  df = df[df$remove == FALSE, ]
  
  dfstat = df[, c('gene.final','rep1','rep2')]
  colnames(dfstat) <- c('gene','rep1','rep2')
  dfstat = dfstat %>% calc_mod_ttest_new() %>% id_significant_proteins()
  
  path.bait.direct = unlist(lapply(strsplit(path,'\\.|(vs)'), function(x) x[2]))
  path.bait = translate_path_bait(path.bait.direct, path_bait_table)
  path.bait.inweb <- gsub('PLPP3', 'PPAP2B', gsub('JCAD', 'KIAA1462',path.bait))
  path.vs = unlist(lapply(strsplit(path,'\\.|(vs)'), function(x) x[3]))
  path.vs[path.vs == 'Whitehead'] <- NA
  path.cell = unlist(lapply(strsplit(path, '\\.'), function(x) x[3]))
  path.date = unlist(lapply(strsplit(path, '\\.'), function(x) x[5]))
  path.facility = ifelse(grepl('Broad', path), 'Broad', 'Whitehead')
  
  print(paste0(path.bait,' (', path.bait.direct, ')...'))
  
  # quality of IP
  ip.repcor <- cor(dfstat$rep1, dfstat$rep2)
  ip.flag.tag <- grepl(tolower(path), 'flag')
  ip.flag.tag.found <- any(!is.na(df$flag.tag))
  #ip.mutant <- any(grepl, df$mutant)
  #ip.wildtype <- any(grepl, df$wildtype)
  ip.median.peptides <- median(df$unique.peptides)
  ip.median.unique.peptides <- median(df$unique.peptides)
  ip.imputed <- sum(df$imputed)
  
  # check bait and interactors (FDR < 0.1, logFC > 0)
  ip.bait.found <- path.bait %in% dfstat$gene | ip.flag.tag.found
  ip.bait.enriched <- ifelse(path.bait %in% dfstat$gene[dfstat$significant], T, F) 
  ip.significant = sum(dfstat$significant)
  ip.total = length(dfstat$significant)
  inweb.list <- get_inweb_list(standardize(path.bait.inweb))
  
  # FDR < 0.1, logFC > 0
  inweb.known.interactors <- as.character(inweb.list[inweb.list$significant, ]$gene)
  ip.known.interactors <- sum(inweb.known.interactors %in% dfstat$gene)
  ip.roselli.all <- sum(roselli.gene %in% dfstat$gene)
  ip.roselli.sig <- sum(roselli.gene %in% dfstat$gene[dfstat$significant])
  intersectDf = data.frame(intersectN = F)
  ip.roselli.stats <- suppressWarnings(calc_hyper(dfstat, data.frame(gene = roselli.gene, significant = T), intersectDf))
  ip.inweb.stats = calc_hyper(dfstat, inweb.list, bait = path.bait)
  
  # FDR < 0.1, LogFC < 0
  dfstat_neg = dfstat %>% calc_mod_ttest_new() %>% id_significant_proteins(fdr_cutoff = 0.1, logfc_dir = 'negative')
  ip.significant.neg <- sum(dfstat_neg$significant)
  ip.bait.enriched.neg <- ifelse(path.bait %in% dfstat_neg$gene[dfstat_neg$significant], T, F) 
  ip.known.interactors.neg <- sum(inweb.known.interactors %in% dfstat_neg$gene)
  ip.roselli.sig.neg <- sum(roselli.gene %in% dfstat_neg$gene[dfstat_neg$significant])
  ip.roselli.stats.neg <- suppressWarnings(calc_hyper(dfstat_neg, data.frame(gene = roselli.gene, significant = T), intersectDf))
  ip.inweb.stats.neg = calc_hyper(dfstat_neg, inweb.list, bait = path.bait)
  
  # FDR < 0.1
  dfstat_both = dfstat %>% calc_mod_ttest_new() %>% id_significant_proteins(fdr_cutoff = 0.1, logfc_dir = 'both')
  ip.significant.both <- sum(dfstat_both$significant)
  ip.bait.enriched.both <- ifelse(path.bait %in% dfstat_both$gene[dfstat_both$significant], T, F) 
  ip.known.interactors.both <- sum(inweb.known.interactors %in% dfstat_both$gene)
  ip.roselli.sig.both <- sum(roselli.gene %in% dfstat_both$gene[dfstat_both$significant])
  ip.inweb.stats.both = calc_hyper(dfstat_both, inweb.list, bait = path.bait)
  ip.roselli.stats.both <- suppressWarnings(calc_hyper(dfstat_both, data.frame(gene = roselli.gene, significant = T), intersectDf))
  

  rowsummary <- data.frame(date = path.date, path = path, 
                            facility = path.facility, 
                            bait = path.bait, 
                            control = path.vs, 
                            cell = path.cell, 
                            imputed = ip.imputed,
                            correlation = ip.repcor, 
                            flag.tag.found = ip.flag.tag.found, 
                            rows.kept.after.qc = paste(rows.qc.discarded, '/', rows.qc.total),
                            peptides.median = ip.median.peptides,
                            unique.peptides.median = ip.median.unique.peptides, 
                            
                            bait.found = ip.bait.found,
                            
                            # FDR < 0.1, logFC > 0
                            bait.enriched.pos = ip.bait.enriched,
                            proteins.significant.pos = paste(ip.significant, '/', ip.total), 
                            proteins.in.inweb.pos = paste(ip.known.interactors, '/', length(inweb.known.interactors), ' p-value = ', round(ip.inweb.stats$statistics$pvalue, 4)),
                            proteins.roselli.pos = paste(ip.roselli.sig, '/', ip.roselli.all, ' p-value = ', round(ip.roselli.stats$statistics$pvalue, 4)),
                            
                            # FDR < 0.1
                            bait.enriched = ip.bait.enriched.both,
                            proteins.significant = paste(ip.significant.both, '/', ip.total), 
                            proteins.in.inweb = paste(ip.known.interactors.both, '/',  length(inweb.known.interactors), ' p-value = ', round(ip.inweb.stats.both$statistics$pvalue, 4)),
                            proteins.roselli = paste(ip.roselli.sig.both, '/', ip.roselli.all, ' p-value = ', round(ip.roselli.stats.both$statistics$pvalue, 4)),
                            
                            # FDR < 0.1, logFC < 0
                            bait.enriched.neg = ip.bait.enriched.neg,
                            proteins.significant.neg = paste(ip.significant.neg, '/', ip.total), 
                            proteins.in.inweb.neg = paste(ip.known.interactors.neg, '/',  length(inweb.known.interactors), ' p-value = ', round(ip.inweb.stats.neg$statistics$pvalue, 4)),
                            proteins.roselli.neg = paste(ip.roselli.sig.neg, '/', ip.roselli.all, ' p-value = ', round(ip.roselli.stats.neg$statistics$pvalue, 4))
                           
                           
                           )

  summary[[path]] <- rowsummary
  
  df_new <- df[,c('gene.final', 'accession','unique.peptides','imputed', 'chrom', 'bp.start', 'bp.end')]
  colnames(df_new) <- c('gene.delete', 'accession_number','unique.peptides','imputed', 'chrom', 'bp.start', 'bp.end')
  df_merge = cbind(dfstat, df_new)
  
  if (!all(dfstat$gene == df_new$gene.delete)) stop('gene rows does not match!')
  
  #df_merge = merge(dfstat, df_new, by = 'gene')
  
  df_merge$gene.delete <- NULL
  
  #if (nrow(df_merge_old) != nrow(df_merge)) stop('data.frames not equal!')
  
  #print(colnames(df_merge))
  
  df_merge <- df_merge[with(df_merge, order(-logFC, FDR)),] 
  # write data 
  
  ### SUBSET COLUMNS (NB: overwrites previous files!)
  write.table(df_merge, file = newpath, sep = '\t', quote = F, row.names = F)
  
  # make volcano plots
  bait_synonyms <- unique(c(path.bait, path.bait.direct, path.bait.direct))
  if (any(bait_synonyms %in% df_merge$gene)) print('bait present!!')
  
  title = paste0(path.bait,' vs ',path.vs, ' in ', path.cell)
  outpath = paste0('~/Projects/03_MICOM/derived/',targetrun,'/',basename(tools::file_path_sans_ext(path)), '.pdf')
  
  # stats
  prot <- rowsummary$proteins.significant.pos
  inweb <- rowsummary$proteins.in.inweb.pos
  roselli <- rowsummary$proteins.roselli.pos
  imputed <- ip.imputed
  
  subtitle <- paste0('Path: ', basename(path), '\n',
                     'Run: ', writerun, '\n',
                     'r: ', ip.repcor, '\n',
                     'Flag Tag: ', ip.flag.tag,'\n',
                     'Enriched:  ', prot, '\n',
                     'Imputed:  ', imputed, '\n',
                     'Roselli:  ', roselli, '\n',
                     'InWeb:  ', inweb, '\n')
  
  
  
  pdf(outpath, width = 8, height = 10)
  # plotting all three options
  volcano = df_merge %>% 
    plot_volcano_basic() %>%
    plot_overlay(as.bait(bait_synonyms))
  
  # Plot different volcano plots
  volcano = volcano + ggtitle(title, subtitle)
  print(volcano)

  volcano_inweb = volcano %>% plot_overlay(list(inweb = inweb.list))
  print(volcano_inweb)
  
  volcano_roselli = volcano %>% plot_overlay(as.goi(roselli.gene, dataset = 'Roselli Genes', col_significant = 'brown'))
  print(volcano_roselli)
  
  # plot scatter plots
  scatter = df_merge  %>% 
    plot_scatter_basic() %>%
    plot_overlay(as.bait(bait_synonyms)) + ggtitle(title)
  
  print(scatter)
  
  graphics.off()
  
  
}


# for renaming data.frame
newcolnames =  c('Analysis Date', 'Data path', 'Facility', 'Bait', 'Control', 
    'Cell line', 'Imputed protein intensities (n)', 'Correlation (r)', 'Flag tag (T/F)', 'Rows after QC (n)', 
    'Peptides found in LC-MS/MS (median)', 'Unique petides found in LC-MS/MS (median)', 'Bait found (T/F)', 
    'Bait enriched (LogFC > 0, FDR < 0.1)','Proteins enriched (LogFC > 0, FDR < 0.1)', 'InWeb interactors enriched (logFC > 0, FDR < 0.1)', 'Roselli genes enriched (LogFC > 0, FDR < 0.1)', 
    'Bait enriched (FDR < 0.1)','Proteins enriched (FDR < 0.1)', 'InWeb interactors enriched (FDR < 0.1)', 'Roselli genes enriched (FDR < 0.1)', 
    'Bait enriched (LogFC < 0, FDR < 0.1)','Proteins enriched (LogFC < 0, FDR < 0.1)', 'InWeb interactors enriched (logFC < 0, FDR < 0.1)', 'Roselli genes enriched (LogFC < 0, FDR < 0.1)' )


main <- do.call(rbind, summary)
main <- main[with(main, order(facility, bait, control, cell, imputed, correlation)),]

# setup tiers
tier1 <- main$correlation >= 0.6 & main$bait.enriched.pos
tier2 <- main$correlation >= 0.4 & main$bait.enriched.pos
tier3 <- main$correlation >= 0.2 & main$bait.enriched.pos
tier4 <- main$correlation >= 0 & main$bait.found
main$tier <- ifelse(tier1, 1, ifelse(tier2, 2, ifelse(tier3, 3, ifelse(tier4, 4, 5))))


prefix


main1 <- main
colnames(main1) <- newcolnames
write.table(main1, paste0('~/Desktop/',prefix,'_summary.tsv'), sep = '\t', quote = F, row.names = F)


main2 <- main[with(main, order(tier, cell, facility, bait, correlation)),]
colnames(main2) <- c(newcolnames, 'tier')
write.table(main2, paste0('~/Desktop/',prefix,'_tier_summary.tsv'),sep = '\t', quote = F, row.names = F)




lst.genome <- list()
# get genomic mapping
for (path in paths){
  
  df = read.csv(path, sep = '\t', stringsAsFactors = F) 
  df = df[df$remove == FALSE, ]
  mydf <- df[,c('gene.final', 'chrom', 'bp.start', 'bp.end')]
  colnames(mydf)[1] <- 'gene'
  lst.genome[[path]] <- mydf
  
  
}

genome <- as.data.frame(do.call(rbind, lst.genome))
rownames(genome) <- NULL
genome <- genome[!duplicated(genome), ]
genome <- genome[with(genome, order(chrom, gene)), ]
write.csv(genome, paste0('~/Desktop/',prefix,'_gene_mapping.csv'), quote = F, row.names = F)

# how many of our baits are associated with ribosomal proteins?

ribo <- lapply(baits, function(b) {
  b <- gsub('PLPP3', 'PPAP2B', gsub('JCAD', 'KIAA1462', b))
  proteins <- get_inweb_list(b)
  proteins <- proteins$gene[proteins$significant]
  return(paste(sum(grepl('^RPL', proteins)),'/',length(proteins)))
})
names(ribo) <- baits








