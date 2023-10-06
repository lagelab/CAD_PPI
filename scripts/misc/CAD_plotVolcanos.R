# micom scripts to generate figure 1 

library(genoppi)
library(data.table)
library(cowplot)

setwd('~/Projects/14_micom_clean/MICOM/')

source('R/read_omim.R')

# ready paths and prefixes
run = 'run056'

# setup paths and run data
files <- list.files('MICOM/data/ms_data_genoppi_input/run056/', full.names = T)


experiments <- read.csv('derived/micom_run056_experiments.csv')
experiments$path <- gsub('data/genoppi_input/','data/ms_data_genoppi_input/',experiments$path)
experiments$bait_mod <- gsub('KIAA1462','JCAD',experiments$bait)


gtex_cat <- read.csv('~/Projects/15_genesets/genesets/data/gtex/GTEX.tstat.categories.genoppi.csv')
cv_tissue = gtex_cat$Tissue.genoppi[gtex_cat$Tissue.category.for.display == 'Cardiovascular']


## InWeb
p1 <- lapply(1:nrow(experiments), function(i){
  
  # gather stats
  df = fread(experiments$path[i])
  bait <- experiments$bait[i]
  bait_db <-  experiments$bait_mod[i]
  db <- get_inweb_list(bait)
  stats <- calc_hyper(df, db, bait = bait, intersectDf = data.frame(intersectN = F))
  
  # generate plot
  plt <- plot_volcano_basic(df) %>% 
    plot_overlay(as.bait(bait_db), label_size = 3) %>% 
    plot_overlay(list(inweb = db)) %>%
    theme_volcano()
  
  # plot details
  title = paste(experiments$bait_mod[i], experiments$cell[i])
  subtitle = paste(experiments$How[i], experiments$Facility[i], '\nInWeb P-value:',round(stats$statistics$pvalue, 4))
  
  
  plt <- plt +  
    ggtitle(title, subtitle) +
    theme(legend.position = "none")
  return(plt)
})
p1[[21]] <- ggplot() + annotate("text", x = 4, y = 25, size=5, label = 'InWeb') + theme_void()

## IRefIndex
p2 <- lapply(1:nrow(experiments), function(i){
  
  # gather stats
  df = fread(experiments$path[i])
  bait <- experiments$bait[i]
  bait_db <-  experiments$bait_mod[i]
  db <- get_irefindex_list(bait)
  if (!is.null(db)){
    stats <- calc_hyper(df, db, bait = bait, intersectDf = data.frame(intersectN = F))
    pvalue = stats$statistics$pvalue
  } else {
    pvalue = 1
  }
  
  # generate plot
  plt <- plot_volcano_basic(df) %>% 
    plot_overlay(as.bait(bait_db), label_size = 3) %>% 
    plot_overlay(list(irefindex = db)) %>%
    theme_volcano()
  
  # plot details

  title = paste(experiments$bait_mod[i], experiments$cell[i])
  subtitle = paste(experiments$How[i], experiments$Facility[i], '\nIRefIndex P-value:',round(pvalue, 4))
  
  
  plt <- plt +  
    ggtitle(title, subtitle) +
    theme(legend.position = "none")
  return(plt)
})
p2[[21]] <- ggplot() + annotate("text", x = 4, y = 25, size=5, label = 'IRefIndex') + theme_void()


## Tissue-specific genes
p3 <- lapply(1:nrow(experiments), function(i){
  
  # gather stats
  df = fread(experiments$path[i])
  bait <- experiments$bait[i]
  bait_db <-  experiments$bait_mod[i]
  db <- data.frame(get_tissue_lists(cv_tissue[2], gtex_rna), col_significant = 'green')
  
  stats <- calc_hyper(df, db, bait = bait, intersectDf = data.frame(intersectN = F))
  
  # generate plot
  plt <- plot_volcano_basic(df) %>% 
    plot_overlay(as.bait(bait), label_size = 3) %>% 
    plot_overlay(list(gtex = db)) %>%
    theme_volcano()
  
  # plot details
  title = paste(experiments$bait_mod[i], experiments$cell[i])
  subtitle = paste(experiments$How[i], experiments$Facility[i], '\nP-value:',round(stats$statistics$pvalue, 4))
  
  
  plt <- plt +  
    ggtitle(title, subtitle) +
    theme(legend.position = "none")
  return(plt)
})
p3[[21]] <- ggplot() + annotate("text", x = 4, y = 25, size=5, label = paste('GTEx ',cv_tissue[2],'\n enrichment')) + theme_void()


## OMIM
ompath <- '~/Projects/15_genesets/genesets/data/omim/genemap2.txt'
re2 = ''
omim <- list(
  stroke = read_omim('stroke', re2, ompath),
  myocardial = read_omim('myocardial', re2, ompath),
  hypertension = read_omim('hypertension', re2, ompath),
  diabetes = read_omim('diabetes', re2, ompath),
  cardiomyopathy = read_omim('Cardiomyopathy', re2, ompath),
  vascular = read_omim('Vascular', re2, ompath)
)


p4<- lapply(1:nrow(experiments), function(i){
  
  # gather stats
  df = fread(experiments$path[i])
  bait <- experiments$bait[i]
  bait_db <-  experiments$bait_mod[i]
  db <- omim$cardiomyopathy
  db$col_significant <- 'blue'
  
  stats <- calc_hyper(df, db, bait = bait, intersectDf = data.frame(intersectN = F))
  
  # generate plot
  plt <- plot_volcano_basic(df) %>% 
    plot_overlay(as.bait(bait), label_size = 3) %>% 
    plot_overlay(list(gtex = db)) %>%
    theme_volcano()
  
  # plot details
  title = paste(experiments$bait_mod[i], experiments$cell[i])
  subtitle = paste(experiments$How[i], experiments$Facility[i], '\n P-value:',round(stats$statistics$pvalue, 4))
  
  
  plt <- plt +  
    ggtitle(title, subtitle) +
    theme(legend.position = "none")
  return(plt)
})
p4[[21]] <- ggplot() + annotate("text", x = 4, y = 25, size=5, label = 'OMIM Cardiomypathy\nenrichment') + theme_void()

pdf('derived/plots/210415_run056_volcanos_inweb.pdf', width = 17, height = 12)
plot_grid(plotlist = p1, nrow = 3, ncol = 7, labels = LETTERS[1:20])
plot_grid(plotlist = p2, nrow = 3, ncol = 7, labels = LETTERS[1:20])
plot_grid(plotlist = p3, nrow = 3, ncol = 7, labels = LETTERS[1:20])
plot_grid(plotlist = p4, nrow = 3, ncol = 7, labels = LETTERS[1:20])
graphics.off()
