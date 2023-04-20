# what genesets are enriched in the WT versus MT.
setwd('~/Projects/14_micom_clean/MICOM/')
library(ggplot2)
library(genoppi)
library(cowplot)
library(writexl)
library(data.table)
source('R/make_cell_heatmap.R')
source('R/read_omim.R')
run = 'run056'

# helpers
index <- function(x, i) x[i]
background_genes <- unique(c(inweb_table$Gene1, inweb_table$Gene2))

# setup paths and run data
gtex_cat <- fread('~/Projects/15_genesets/genesets/data/gtex/GTEX.tstat.categories.genoppi.csv')
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

# manually import networks
read_network <- function(d){return(data.frame(gene = d$GeneName, significant = d$Interactor))}

EC <- read_network(fread('data/interactor_lists/run056/COMBINED_EC.InteractorTable.txt'))
SMC <- read_network(fread('data/interactor_lists/run056/COMBINED_SMC.InteractorTable.txt'))
EC_SMC <- read_network(fread('data/interactor_lists/run056/COMBINED_EC-SMC.InteractorTable.txt'))

# EC/SMC networks
data[['EC']] = EC
data[['SMC']] = SMC
data[['EC_SMC']] = EC_SMC

# Tissue Specific networks
data[['EC_only']] =  EC[! EC$gene %in% SMC$gene,]
data[['SMC_only']] =  SMC[! SMC$gene %in% EC$gene,]
data[['EC_SMC_intersect']] = EC_SMC[EC_SMC$gene %in% EC$gene & EC_SMC$gene %in% SMC$gene,]


# counts
#sum(data$EC$significant)
#sum(data$EC$significant)
sum(data$EC$gene[data$EC$significant] %in% data$SMC$gene[data$SMC$significant])
sum(data$EC$significant)

#sum(data$SMC$significant)
#sum(data$SMC$significant)
sum(data$SMC$gene[data$SMC$significant] %in% data$EC$gene[data$EC$significant])
sum(data$SMC$significant)



hyper <- calc_hyper(data$EC, data$SMC)
# 
pdf('derived/plots/210422_ec_smc_venn_diagram.pdf',width = 6, height = 6)
title = paste('EC-SMC overlap\nP-value = 1.4854e-24')
diagram <- list(EC = data$EC$gene[data$EC$significant], 
                SMC = data$SMC$gene[data$SMC$significant])
venn <- draw_genoppi_venn(diagram, colors = c('red','blue'), main =  title)
grid::grid.draw(venn)
graphics.off()





