# what genesets are enriched in the WT versus MT.
setwd('~/Projects/14_micom_clean/MICOM/')
library(ggplot2)
library(genoppi)
library(reshape2)
library(writexl)
library(data.table)
source('R/make_cell_heatmap.R')
run = 'run056'

# helpers
index <- function(x, i) x[i]
background_genes <- unique(c(inweb_table$Gene1, inweb_table$Gene2))

# haarst genes
haarst <- read_excel('derived/genelists/27AUG20_haarst2017_ensembl_anneal_mapping_tags.xlsx')
haarst_genes <- unique(na.omit(haarst$gene_name))

# setup paths and run data
files <- list.files('MICOM/data/ms_data_genoppi_input/run056/', full.names = T)
experiments <- read.csv('derived/micom_run056_experiments.csv')
experiments$path <- gsub('data/genoppi_input/','data/ms_data_genoppi_input/',experiments$path)
experiments$submission <- lapply(strsplit(basename(experiments$path), split = '\\.'), index, 2)
experiments$id <- paste(experiments$bait,
                        experiments$cell,
                        ifelse(experiments$How == 'Overexpression', 'OE','EN'), 
                        ifelse(experiments$Facility == 'Broad', 'BR','WH'),
                        toupper(experiments$submission), sep = '_')

# load IP-MS data 
data <- lapply(experiments$path, function(x) fread(x))
names(data) <- experiments$id
baits <- unique(experiments$bait)
cells <- unique(experiments$cell)

# get data where we have one SMC and one EC
cell_bait = as.matrix(table(experiments$bait, experiments$cell))
bool <- apply(cell_bait, 1, all)
baits_ok <- rownames(cell_bait)[bool]
experiments <- experiments[experiments$bait %in% baits_ok,]

# load each bait in each cell type and compare them
lst <- list()
for (bait in baits_ok){
  
  # regex the right IDs
  bool_bait <- grepl(bait, experiments$bait)
  row_smc <- experiments[bool_bait  & experiments$cell == 'SMC',][1,]
  row_ec <- experiments[bool_bait  & experiments$cell == 'EC',][1,]
  
  
  # we only compare OE to OE, and Endogenous to Endogenous
  if (row_smc$How == row_ec$How){
    cell_smc <- row_smc$id #experiments$id[bool_bait  & experiments$cell == 'SMC'][1]
    cell_ec <- row_ec$id #experiments$id[bool_bait  & experiments$cell == 'EC'][1]
      
      # load the files.
      data_smc <- data[names(data) %in% cell_smc][[1]]
    data_ec <-  data[names(data) %in% cell_ec][[1]]
    
    # compare the two files
    colnames(data_smc) <- paste0('SMC.',colnames(data_smc))
    colnames(data_ec) <- paste0('EC.',colnames(data_ec))
    data_mrg <- merge(data_ec, data_smc, by.x = 'EC.gene', by.y = 'SMC.gene')
    data_mrg$logd <- data_mrg$SMC.logFC - data_mrg$EC.logFC # log scale so minus
    
    # get known interactors
    inweb <- get_inweb_list(bait)
    sum(inweb$significant)
    data_mrg$inweb <- data_mrg$EC.gene %in% inweb$gene[inweb$significant]
    
    # setup labels
    data_mrg$how <- row_smc$How
    data_mrg$bait <- bait
    data_mrg$bait_label <- apply(data_mrg[,c('bait','how')],1, paste, collapse = '\n')
    
    
    data_mrg$significant <- ifelse(data_mrg$EC.significant & data_mrg$SMC.significant, '***',
                                   ifelse(data_mrg$EC.significant | data_mrg$SMC.significant, '**', '*'))
  
    lst[[bait]] <- data_mrg
  }
}


# get the data
res <- do.call(rbind, lst)
restable <- res

# remvoe some data
res <- res[! res$bait %in% 'HDAC9', ]
differential <- res[abs(res$logd) > 0.7,]
res <- res[res$EC.gene %in% differential$EC.gene,]

# convert to wide for clustering
df_wide <- reshape(res[,c('EC.gene','bait','logd')], 
                   direction = "wide",
                   v.names = "logd",
                   idvar = "EC.gene",
                   timevar = "bait")

df_wide[is.na(df_wide)] <- 0
tree <- hclust(dist(df_wide[,-1]), method = 'complete')
df_wide$cluster_complete_order <- tree$order
levels <- df_wide$EC.gene[tree$order]

# set order to clustered order
plot_df <- res
plot_df$gene <- factor(plot_df$EC.gene, levels = levels)
plot_df$grid <- 'SMC / EC'

# Various plotting stats

# 
plot_df$label <- ifelse(plot_df$inweb, 'InWeb','')
ylabel_cols <- ifelse(levels %in% baits, "red", 
                      ifelse(levels %in% haarst_genes, "blue", "black"))

title = 'Differential interactome for IP experiments (Complete Clustering)'
subtitle = 'Interactor inclusion criteria: at least one IP has abs(log2(SMC/EC)) > 0.7)'

# plot result
graphics.off()
ggplot(plot_df, aes(x = bait_label, y = gene, fill = logd, label = significant)) +
  scale_fill_gradient2(low = 'yellow', mid = 'white', high = 'red', name = 'Log2(SMC/EC)') +
  geom_tile() + 
  geom_text(size = 2) + 
  ggtitle(title, subtitle) +
  theme_bw() + xlab('IP-MS Bait') + ylab('Prey') +
  facet_wrap(~grid) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour = ylabel_cols, size = 5),
        axis.text.x = element_text(angle = 0, size = 8))
ggsave('derived/plots/210416_cell_specificity_heatmap.pdf', width = 7, height = 24)






###### ONLY SELECETED LOCI

# get the data
res <- do.call(rbind, lst)
res <- res[! res$bait %in% 'HDAC9', ]
differential <- res[res$EC.gene %in% haarst_genes | res$EC.gene %in% baits_ok,]
res <- res[res$EC.gene %in% differential$EC.gene,]

# convert to wide for clustering
df_wide <- reshape(res[,c('EC.gene','bait','logd')], 
                   direction = "wide",
                   v.names = "logd",
                   idvar = "EC.gene",
                   timevar = "bait")

df_wide[is.na(df_wide)] <- 0
tree <- hclust(dist(df_wide[,-1]), method = 'complete')
df_wide$cluster_complete_order <- tree$order
levels <- df_wide$EC.gene[tree$order]

# set order to clustered order
plot_df <- res
plot_df$gene <- factor(plot_df$EC.gene, levels = levels)
plot_df$grid <- 'SMC / EC'

# Various plotting stats

# 
plot_df$label <- ifelse(plot_df$inweb, 'InWeb','')
ylabel_cols <- ifelse(levels %in% baits, "red", 
                      ifelse(levels %in% haarst_genes, "blue", "black"))

title = 'Differential interactome for IP experiments (Complete Clustering)'
subtitle = 'Interactors that are in Haarst 2017 tagged loci and baits'

# plot result
graphics.off()
ggplot(plot_df, aes(x = bait_label, y = gene, fill = logd, label = significant)) +
  scale_fill_gradient2(low = 'yellow', mid = 'white', high = 'red', name = 'Log2(SMC/EC)') +
  geom_tile() + 
  geom_text(size = 4) + 
  ggtitle(title, subtitle) +
  theme_bw() + xlab('IP-MS Bait') + ylab('Prey') +
  facet_wrap(~grid) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(colour = ylabel_cols),
        axis.text.x = element_text(angle = 0))
ggsave('derived/plots/210416_cell_specificity_selected_heatmap.pdf', width = 7, height = 10)





##  get res table
d <- restable
d <- d[,c(29, 1:6,14:18,27,19,7,26)]
d$logd_gt_one <- ifelse(abs(d$logd) > 1, 'T','F')
d <- d[rev(order(abs(d$logd))), ]
d <- d[!duplicated(d),]
write.csv(d, '~/Desktop/210418_supplementary_table5.csv', quote = T)






