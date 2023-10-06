
library(readxl)
devtools::load_all('')
source('R/read_excel_sheets.R')
source('R/strsplit.newline.R')
source('R/ggbarplot.R')

# get gtex
df_gtex <- read_excel_sheets('derived/210422_gtex_rna_global_network_table.xlsx')
gtex_cat <- read.csv('~/Projects/15_genesets/genesets/data/gtex/GTEX.tstat.categories.genoppi.csv')

outfile = 'derived/plots/210423_gtex_rna_network_table_ordered_by_pvalue.pdf'
# gtex plots
pdf(outfile, width = 5, height = 6)
for (dname in names(df_gtex)){
  
  df <- df_gtex[[dname]]
  #df <- merge(cur_df, gtex_cat, by.x = 'dataset', by.y = 'Tissue.genoppi')
  
  # subset by p-value and global enrichment
  df <- df[df$enrichment == 'global',]
  df <- df[df$pvalue != 1,]
  
  # order columns
  #df$Tissue <- factor(df$Tissue, levels = unique(df$Tissue[order(df$Tissue.category.for.display)]))
  df$Tissue <- factor(df$Tissue, levels = unique(df$Tissue[order(df$Tissue.category.for.display, df$pvalue)]))
  df$enrichment <- factor(df$enrichment, levels = c('global','conditional'))
  
  # log-pvalue
  df$logpvalue <- -log10(df$pvalue)
  
  bonf <- 0.05  / length(unique(df$Tissue))
  
  # set significance thresholds
  df$significance <- 'Not significant'
  df$significance[df$pvalue < 0.05] <- 'Nominally significant'
  df$significance[df$pvalue < bonf] <- 'Bonferroni Corrected'
  df$significance <- as.factor(df$significance)
  
  # generate color palette
  #colors = c("firebrick2","tomato1","black")
  #names(colors) <- c('FDR < 0.1','Nominally significant', 'Not significant')
  #color_scale <- scale_fill_manual(name = "significance", values = colors)

  # generate plot
  plt <- ggplot(df, aes(y = Tissue, x = logpvalue, fill = Tissue.category.for.display)) + 
    geom_bar(stat = 'identity', alpha = 0.9) +
    geom_vline(xintercept = -log10(0.05), color = "black", linetype = 'dashed') + 
    geom_vline(xintercept = -log10(bonf), color = "black", linetype = 'dashed') + 
    xlab('-log10(P-value)') + ylab('GTEx tissue') +
    ggtitle(paste(dname, '(GTEx - RNA data)')) + 
    labs(fill='Tissue Category') +
    theme_bw() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

  plot(plt)
}
graphics.off()





##### ALL GTEX TISSUE #####

# individual networks combined
df <- do.call(rbind,lapply(names(df_gtex)[c(1:6)], function(dname){
  df <- df_gtex[[dname]]
  #df <- merge(cur_df, gtex_cat, by.x = 'dataset', by.y = 'Tissue.genoppi')
  df$name <- dname
  return(df)
}))

# order columns
#df$Tissue <- factor(df$Tissue, levels = unique(df$Tissue[order(df$Tissue.category.for.display)]))
df <- df[df$enrichment == 'global',]
df$list_name[df$list_name == 'EC_SMC_union'] <- 'Union'
df$list_name[df$list_name == 'EC_SMC_intersect'] <- 'Intersect'
df$list_name[df$list_name == 'EC_only'] <- 'EC only'
df$list_name[df$list_name == 'SMC_only'] <- 'SMC only'


df$list_name<- factor(df$list_name, levels = c('EC','Union','SMC','EC only','Intersect', 'SMC only'))
df$Tissue <- factor(df$Tissue, levels = unique(df$Tissue[order(df$Tissue.category.for.display, df$pvalue)]))

# log-pvalue
df$logpvalue <- -log10(df$pvalue)
bonf <- 0.05  / length(unique(df$Tissue))

# generate plot
outfile = 'derived/plots/210424_gtex_rna_across_all.pdf'
pdf(outfile, width = 6, height = 8)
ggplot(df[df$enrichment == 'global',], aes(y = Tissue, x = logpvalue, fill = Tissue.category.for.display)) + 
  geom_bar(stat = 'identity', alpha = 0.9) +
  geom_vline(xintercept = -log10(0.05), color = "black", linetype = 'dashed') + 
  geom_vline(xintercept = -log10(bonf), color = "black", linetype = 'dashed') + 
  xlab('-log10(P-value)') + ylab('GTEx tissue') +
  ggtitle('Integrated PPI Networks Tissue Global Enrichment (GTEx - RNA data)') + 
  labs(fill='Tissue Category') +
  theme_bw() +
  theme(legend.position='bottom',plot.title=element_text(size=9),
        axis.title=element_text(size=9),axis.text=element_text(size=7),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_wrap(~list_name)
graphics.off()




##### CARDIOVASCULAR GTEX TISSUE #####

# selected tissues 
ndf <- df[df$enrichment == 'global' & df$Tissue.category.for.display %in% 'Cardiovascular',]
ndf$experiment <- ndf$list_name
ndf$experiment <- factor(ndf$experiment, levels = c('EC','EC_SMC_union','SMC','EC_only','EC_SMC_intersect', 'SMC_only'))
ndf$list_name <- ndf$Tissue


outfile = 'derived/plots/210425_gtex_rna_Cardiovascular.pdf'
pdf(outfile, width = 1.4*2, height = 2)
for (selection in unique(ndf$experiment)){
  
  plt <- ggbarplot(ndf[ndf$experiment %in% selection,], bonf = (0.05 / 53)) +  
    ggtitle(gsub('_',' ',selection)) + theme_bw() + 
    theme(legend.position='none',plot.title=element_text(size=9),
          axis.title=element_text(size=9),axis.text=element_text(size=7))
  print(plt)
    
}
graphics.off()

#outfile = 'derived/plots/210423_gtex_rna_Cardiovascular.pdf'
#pdf(outfile, width = 8, height = 5)
#ggbarplot(ndf, bonf) + facet_wrap(~experiment)
#ggbarplot(ndf[ndf$experiment %in% c('EC','SMC','EC_SMC'),], bonf,  ylab = 'Cardiovascular tissue') + facet_wrap(~experiment)
#graphics.off()

# individual networks
outfile = 'derived/plots/210423_gtex_rna_indnetworks_Cardiovascular.pdf'
pdf(outfile, width = 5, height = 2)
ggbarplot(ndf[ndf$experiment == 'EC',], bonf, ylab = 'GTEx cardiovascular tissue') + ggtitle('EC Network')
ggbarplot(ndf[ndf$experiment == 'SMC',], bonf, ylab = 'GTEx cardiovascular tissue') + ggtitle('SMC Network')
ggbarplot(ndf[ndf$experiment == 'EC_SMC',], bonf, ylab = 'GTEx Cardiovascular tissue') + ggtitle('EC+SMC Network')
graphics.off()
# individual plots



### Reactome and GO MF enrichment
generate_excel_barplots <- function(path, outfile, height, width, top = 10, label_nchar_split = 30){
  geneset <- read_excel_sheets(path)
  pdf(outfile, width = width, height = height)
  for (g in geneset){
    for (e in c('global','conditional')){
      d <- head(g[g$enrichment == e,], top)
      title <- paste(unique(d$list_name), e, sep = ' - ')
      d$list_name <- unlist(strsplit.newline(d$dataset, label_nchar_split))
      bonf <- 0.05 / length(d$dataset)
      plt <- ggbarplot(d, bonf) +
        ggtitle(title) +
        theme(legend.position='none',plot.title=element_text(size=9),
              axis.title=element_text(size=9),axis.text=element_text(size=7))
      print(plt)
    }
  }
  graphics.off()
}


path <- 'derived/210422_reactome_network_only_table.xlsx'
outfile <- 'derived/210422_reactome_network_only_table.pdf'
generate_excel_barplots(path, outfile, 3.5, 3, label_nchar_split = 30)

path <- 'derived/210422_mf_network_only_table.xlsx'
outfile <- 'derived/210422_mf_network_only_table.pdf'
generate_excel_barplots(path, outfile, 3, 3,  label_nchar_split = 36)

path <- 'derived/210422_msigdb_h_network_only_table.xlsx'
outfile <- 'derived/210422_msigdb_h_network_only_table.pdf'
generate_excel_barplots(path, outfile, 5, 5, label_nchar_split = 35)

path <- 'derived/210422_bp_network_only_table.xlsx'
outfile <- 'derived/210422_bp_network_only_table.pdf'
generate_excel_barplots(path, outfile, 3, 3, label_nchar_split = 35)

path <- 'derived/210422_cc_network_only_table.xlsx'
outfile <- 'derived/210422_cc_network_only_table.pdf'
generate_excel_barplots(path, outfile, 3, 3, label_nchar_split = 35)

