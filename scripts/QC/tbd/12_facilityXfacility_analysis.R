setwd('~/Projects/03_MICOM/')

#comparisons:
# (1) EC (7 IPs) vs. SMC (11 IPs)
# (2) Broad (3 IPs) vs. Whitehead (15 IPs)
# (3) overexpression (5 IPs) vs. endogenous (13 IPs)
# metrics to show in boxplots:
# (1) replicate correlation
# (2) # detected proteins
# (3) # significant proteins
# (4) % significant proteins (significant/detected)
# (5) % ribosomal proteins (ribosomal/detected)
# (6) % significant ribosomal proteins (ribosomal/significant)
# (7) -log10 InWeb p-value


### GET DATA

# get summary data
summary.data <- read.csv('201120_run056_micom_tier_summary.tsv', sep = '\t') #summary.data <- read.csv('~/Projects/03_MICOM/21JUL20_run019_micom_summary.csv')
summary.data <- summary.data[,c('Data.path', 'Correlation..r.', 'Peptides.found.in.LC.MS.MS..median.', 'Unique.petides.found.in.LC.MS.MS..median.', 
                                'Proteins.enriched..LogFC...0..FDR...0.1.', 'InWeb.interactors.enriched..logFC...0..FDR...0.1.', 'Roselli.genes.enriched..LogFC...0..FDR...0.1.')]
colnames(summary.data) <- c('Path','r', 'peptides', 'unique.peptides', 'proteins.enriched','inweb','roselli')

# get qc data, i.e. the  IPs we want to investigate
qc.data <- data.frame(Path=read.csv('run056_filtered_tier1_paths.tsv', sep = '\t')$x)  #qc.data <- read.csv('~/Desktop/29JUL20_MICOM_qced_details_FM.csv')
qc.data$Path <- gsub('run056', 'run055', qc.data$Path)
sum(summary.data$Path %in% qc.data$Path) # expect 20..
summary.data <- summary.data[summary.data$Path %in% qc.data$Path,]
nrow(summary.data)
# lysate input (mg) should be its own column
#lysate <- strsplit(qc.data$Lysate.Input, ' x ')
#qc.data$lysate.input.mg <- as.numeric(unlist(gsub('mg', '', lapply(lysate, function(x) x[1]))))
#qc.data$lysate.plex <- unlist(lapply(lysate, function(x) x[2]))
#qc.data$Lysate.Input <- NULL

# get ribosomal proteins
ribosomal.proteins <- lapply(summary.data$Path, function(p){
  p = gsub('run055', 'run056', p)
  df <- read.csv(p, sep = '\t')
  df$ribosomal <- grepl('^RPL', df$gene) | grepl('^RPS', df$gene)
  ribosomal.total <- sum(df$ribosomal)
  ribosomal.sig <- sum(df$ribosomal & df$significant)
  return(data.frame(Path = p, ribosomal.significant = ribosomal.sig, ribosomal.detected = ribosomal.total, ribosomal.significant.pct = ribosomal.sig / sum(df$significant)))
})
ribosomal.data = do.call(rbind,ribosomal.proteins)

# extract info from summary stats
inweb = strsplit(summary.data$inweb, ' ')
roselli = strsplit(summary.data$roselli, ' ')
proteins = strsplit(summary.data$proteins.enriched, ' ')

# proteins detected / significant
summary.data$proteins.significant <- as.numeric(unlist(lapply(proteins, function(x) x[1] )))
summary.data$proteins.detected <- as.numeric(unlist(lapply(proteins, function(x) x[3] )))
summary.data$proteins.significant.pct <- summary.data$proteins.significant / summary.data$proteins.detected
summary.data$proteins.enriched <- NULL

# inweb detected / significant / pvalue
summary.data$inweb.significant <- as.numeric(unlist(lapply(inweb, function(x) x [1])))
summary.data$inweb.detected <- as.numeric(unlist(lapply(inweb, function(x) x [3])))
summary.data$inweb.pvalue.nlog10 <- -log10(as.numeric(unlist(lapply(inweb, function(x) x [8]))))
summary.data$inweb.significant.pct <- summary.data$inweb.significant / summary.data$inweb.detected
summary.data$inweb <- NULL

# roselli genes detected / significant / pvalue
summary.data$roselli.significant <- as.numeric(unlist(lapply(roselli, function(x) x [1])))
summary.data$roselli.detected <- as.numeric(unlist(lapply(roselli, function(x) x [3])))
summary.data$roselli.pvalue.nlog10 <- -log10(as.numeric(unlist(lapply(roselli, function(x) x [8]))))
summary.data$roselli.significant.pct <- summary.data$roselli.significant / summary.data$roselli.detected
summary.data$roselli <- NULL


# combine data 
qc.data.corr <- merge(qc.data, summary.data, all.x = T)


#qc.data.corr[,c(1,2, 16)]

# get their tag
qc.data.corr$Tag <- NA
tag.flag <- grepl('3xflag',tolower(qc.data.corr$Path))
tag.endo <- grepl('endoge',tolower(qc.data.corr$Path))
qc.data.corr$Tag[tag.flag] <- 'Overexpression'
qc.data.corr$Tag[tag.endo] <- 'Endogenous'
qc.data.corr$Tag[qc.data.corr$Path == 'data/genoppi_input/run055/Broad.mB.ADAMTS7vsMock.SMC.expanded.20NOV2020.tsv'] <- 'Endogenous'
qc.data.corr$Tag[qc.data.corr$Path == 'data/genoppi_input/run055/Broad.mB.EDN1vsMock.SMC.expanded.20NOV2020.tsv'] <- 'Overexpression'
qc.data.corr$Tag[qc.data.corr$Path == 'data/genoppi_input/run055/Broad.mB.ARHGEF26vsMock.EC.nonSGS.20NOV2020.tsv'] <- 'Overexpression'


# get cell line
qc.data.corr$Cell.line <- unlist(lapply(strsplit(basename(qc.data.corr$Path), split = '\\.|vs'), function(x) x[5]))

# get facility
qc.data.corr$Facility <- unlist(lapply(strsplit(basename(qc.data.corr$Path), split = '\\.|vs'), function(x) x[1]))


#qc.data.corr$Tag[qc.data.corr$Tag == ''] <- 'TBD'
qc.data.corr$Path <- gsub('run055', 'run056', qc.data.corr$Path)
mat <- merge(qc.data.corr, ribosomal.data, by = 'Path')



# deal ambigious labels
#mat$Tag[mat$Tag %in% c('endogenous', "none; endogenous protein")] <- 'Endogenous'
#mat$Tag[mat$Tag %in% c('3xFlag')] <- 'Overexpression'

### PLOT DATA
library(ggplot2)
library(ggsignif)
require(gridExtra)
require(cowplot)




# signif_level
maplevel <- function(p){
  if (p < 0.05){
    return(sprintf("p = %.2g", p))
  } else {
    return('NS.')
  }
}

# comparisons and colors
my.comparisons = list('Cell.line', 'Facility', 'Tag')
my.color.themes <- list('Cell.line' = 'Dark2', 'Facility' = 'Blues', 'Tag' = 'GnBu')


for (comparison in my.comparisons){
  
  pdf(paste0('201128_MICOM_Analysis_Comparison_',comparison,'.pdf'), width = 16, height = 12)
  
  # (1) replicate correlation
  p1 = ggplot(mat, mapping = aes_string(x=comparison, y = 'r', fill = comparison)) +
    ylab('Replicate correlation (r)') + 
    ggtitle('A') + 
    geom_boxplot() +
    geom_jitter(shape = 1) + 
    geom_signif(comparisons = list(unique(mat[[comparison]])),   map_signif_level = maplevel) + #map_signif_level=function(p)sprintf("p = %.2g", p)) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab('') +
    scale_fill_brewer(palette=my.color.themes[[comparison]])
  
  #print(p1)
  
  # (2) # detected proteins
  p2 = ggplot(mat, mapping = aes_string(x=comparison, y = 'proteins.detected',  fill = comparison)) +
    ylab('Proteins detected (n)') + 
    ggtitle('B') + 
    geom_boxplot() +
    geom_jitter(shape = 1) +  
    geom_signif(comparisons = list(unique(mat[[comparison]])), map_signif_level=TRUE) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab('') +
    scale_fill_brewer(palette=my.color.themes[[comparison]])
  
  #print(p2)
  
  # (3) # significant proteins
  p3 = ggplot(mat, mapping = aes_string(x=comparison, y = 'proteins.significant',  fill = comparison)) +
    ylab('Proteins signficant (n)') + 
    ggtitle('C') + 
    geom_boxplot() +
    geom_jitter(shape = 1) +  
    geom_signif(comparisons = list(unique(mat[[comparison]])), map_signif_level=TRUE) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab('') +
    scale_fill_brewer(palette=my.color.themes[[comparison]])

  
  #print(p3)
  
  # (4) % significant proteins (significant/detected)
  p4 = ggplot(mat, mapping = aes_string(x=comparison, y = 'proteins.significant.pct',  fill = comparison)) +
    ylab('Proteins signficant (%)') + 
    ggtitle('D') + 
    geom_boxplot() +
    geom_jitter(shape = 1) +  
    geom_signif(comparisons = list(unique(mat[[comparison]])), map_signif_level=TRUE) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab('') +
    scale_fill_brewer(palette=my.color.themes[[comparison]])

  
  #print(p4)
  
  # (5) % ribosomal proteins (ribosomal/detected)
  p5 = ggplot(mat, mapping = aes_string(x=comparison, y = 'ribosomal.detected', fill = comparison)) +
    ylab('Ribosomal proteins detected (n)') +  
    ggtitle('E') + 
    geom_boxplot() +
    geom_jitter(shape = 1) +  
    geom_signif(comparisons = list(unique(mat[[comparison]])), map_signif_level=TRUE) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab('') +
    scale_fill_brewer(palette=my.color.themes[[comparison]])
  
  #print(p5)
  
  # (6) ribosomal proteins significant
  p6 = ggplot(mat, mapping = aes_string(x=comparison, y = 'ribosomal.significant', fill = comparison)) +
    ylab('Ribosomal proteins significant (n)') + 
    ggtitle('F') +
    geom_boxplot() +
    geom_jitter(shape = 1) +  
    geom_signif(comparisons = list(unique(mat[[comparison]])), map_signif_level=TRUE) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab('') +
    scale_fill_brewer(palette=my.color.themes[[comparison]])
  
  #print(p6)
  
  # (7) % significant ribosomal proteins (ribosomal/significant)
  p7 = ggplot(mat, mapping = aes_string(x=comparison, y = 'ribosomal.significant.pct', fill = comparison)) +
    ylab('Ribosomal proteins significant (%)') + 
    ggtitle('G') +
    geom_boxplot() +
    geom_jitter(shape = 1) +  
    geom_signif(comparisons = list(unique(mat[[comparison]])), map_signif_level=TRUE) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab('') +
    scale_fill_brewer(palette=my.color.themes[[comparison]])
  
  #print(p7)
    
  # (9) Roselli p-value
  p8 = ggplot(mat, mapping = aes_string(x=comparison, y = 'roselli.pvalue.nlog10', fill = comparison)) +
    ylab('-log10(p-value)') + 
    ggtitle('roselli (p-value)') +
    geom_boxplot() +
    geom_jitter(shape = 1) +  
    geom_signif(comparisons = list(unique(mat[[comparison]])), map_signif_level=TRUE) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab('') +
    scale_fill_brewer(palette=my.color.themes[[comparison]])
    
  #print(p8)
  
  # (9) -log10 InWeb p-value
  p9 = ggplot(mat, mapping = aes_string(x=comparison, y = 'ribosomal.significant.pct', fill = comparison)) +
    ylab('InWeb overlap enrichment (-log10 p-value)') + 
    ggtitle('H') +
    geom_boxplot() +
    geom_jitter(shape = 1) +  
    geom_signif(comparisons = list(unique(mat[[comparison]])), map_signif_level=TRUE) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab('') +
    scale_fill_brewer(palette=my.color.themes[[comparison]])
  
  #print(p9)
  p = plot_grid(p1, p2, p3, p4, 
            p5, p6, p7, p9, ncol = 4, nrow = 2)
 
  print(p)
  
  graphics.off()
}



