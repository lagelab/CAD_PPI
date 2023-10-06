# generate heatmaps for geneset enrichments

library(readxl)
source('R/read_excel_sheets.R')

# plot the heatmap
plot_heatmap <- function(df, enrichment  = 'global', threshold = NULL){
  
  # order by p-valie
  df$logpvalue <- -log10(df$pvalue)
  bonf <- 0.05  / length(unique(df$dataset))
  if (!is.null(threshold)) threshold <- threshold else threshold <- bonf
  
  # factorize columns
  df$list_name <- factor(df$list_name, levels = c('EC','SMC','EC_SMC','EC_only','EC_SMC_intersect', 'SMC_only'))
  df$experiment <- df$list_name
  df$list_name <- df$dataset
  
  # subset data 
  cur_df <- df
  cur_df <- cur_df[cur_df$enrichment == enrichment,]
  keep <- unique(cur_df$list_name[cur_df$pvalue < threshold])
  cur_df <- cur_df[cur_df$list_name %in% keep,]
  
  # setup plot
  cur_df$label <- ifelse(cur_df$pvalue < 0.05, ifelse(cur_df$pvalue < bonf, '**','*'),'')
  plt <- ggplot(cur_df, aes(y = reorder(dataset, -log10(pvalue)), x = experiment, label = label, fill = -log10(pvalue))) +
    geom_tile() + geom_text() + theme_bw() +
    scale_fill_gradient(low = 'white', high = 'red') +
    xlab('PPI network') + 
    ylab('GO term') +
    labs('-log10(P-value)')
  return(plt)
  
}

grab_pathways <- function(df, n = 5){
  top <- df[df$enrichment == 'global',]
  pathways <- unique(unlist(lapply(unique(df$list_name), function(x){
    grab <- top[top$list_name %in% x,]
    grab <- head(grab[order(grab$pvalue),], n)
    return(grab$dataset)
  })))
  return(pathways)
}


# load paths
paths <- list.files('derived/', pattern = '210414', full.names = T)
paths <- paths[grepl('(bp)|(cc)|(mf)',paths,)]

## biological process
path <- paths[1]
geneset <- read_excel_sheets(path)
df <- do.call(rbind, geneset[1:6])
selected_pathways <- grab_pathways(df, 5)

# generate plots
plot_heatmap(df[df$dataset %in% selected_pathways,], enrichment = 'global') + 
  ggtitle(paste(basename(path),'global')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.text.y = element_text(angle = 45))
ggsave('derived/plots/210415_go_bp_enrichment_networks_global.pdf', width = 6, height = 6)
#plot_heatmap(df, enrichment = 'conditional', threshold = 0.005) + ggtitle(paste(basename(path),'conditional'))
#ggsave('derived/plots/210414_go_bp_enrichment_networks_conditional.pdf', width = 10, height = 10)

## cellular compartment
path <- paths[2]
geneset <- read_excel_sheets(path)
df <- do.call(rbind, geneset[1:6])
selected_pathways <- grab_pathways(df, 5)

plot_heatmap(df[df$dataset %in% selected_pathways,], enrichment = 'global') + 
  ggtitle(paste(basename(path),'global')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.text.y = element_text(angle = 45))
ggsave('derived/plots/210415_go_cc_enrichment_networks_global.pdf', width = 6, height = 6)
#plot_heatmap(df, enrichment = 'conditional', threshold = 0.005) + ggtitle(paste(basename(path),'conditional'))
#ggsave('derived/plots/210414_go_cc_enrichment_networks_conditional.pdf', width = 8, height = 10)


## molecular function
path <- paths[3]
geneset <- read_excel_sheets(path)
df <- do.call(rbind, geneset[1:6])
selected_pathways <- grab_pathways(df, 5)

plot_heatmap(df[df$dataset %in% selected_pathways,], enrichment = 'global') + 
  ggtitle(paste(basename(path),'global')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.text.y = element_text(angle = 45))
ggsave('derived/plots/210415_go_mf_enrichment_networks_global.pdf', width = 6, height = 6)
#plot_heatmap(df, enrichment = 'conditional', threshold = 0.005) + ggtitle(paste(basename(path),'conditional'))
#ggsave('derived/plots/210414_go_mf_enrichment_networks_conditional.pdf', width = 8, height = 10)

## msigdb H
path <- "derived/210414_msigdb_h_network_table.xlsx" 
geneset <- read_excel_sheets(path)
df <- do.call(rbind, geneset[1:6])
selected_pathways <- grab_pathways(df, 5)

# generate plots
plot_heatmap(df[df$dataset %in% selected_pathways,], enrichment = 'global',  threshold = 0.05) +
  ggtitle(paste(basename(path),'global')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.text.y = element_text(angle = 45))
ggsave('derived/plots/210415_msigdb_h_enrichment_networks_global.pdf',  width = 6, height = 4)
#plot_heatmap(df, enrichment = 'conditional', threshold = 0.05) + ggtitle(paste(basename(path),'conditional'))
#ggsave('derived/plots/210414_msigdb_h_enrichment_networks_conditional.pdf', width = 10, height = 10)



# generate plot
#outfile = 'derived/plots/210414_mf_network.pdf'
#pdf(outfile, width = 8, height = 12)
##for (name in names(geneset)[1:3]){
#  cur_df <- df[df$experiment == name,]
#  # order columns
#  cur_df <- df
#  cur_df <- cur_df[cur_df$enrichment == 'global',]
#  cur_df <- cur_df[cur_df$FDR < 0.05,]
#  plt <- ggbarplot(cur_df, bonf, colors = cols) +
#    ggtitle(paste('GO molecular function global enrichment -', name)) +
#    ylab('GO molecular function category') + 
#    facet_wrap(~experiment)
#  print(plt)
##}
#graphics.off()

