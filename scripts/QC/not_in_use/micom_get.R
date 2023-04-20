

## remove these functions at some point

get.interactors.genes <- function(data, known.interactors){
  
  known.interactors.lst <- as.character(known.interactors[known.interactors$significant == TRUE, ]$gene)
  in.data = data[data$gene %in% known.interactors.lst, ]
  interactors.genes <- NA
  if (nrow(in.data) > 0){
    interactors.ok <- designate(in.data, FDR < 0.1, logFC > 0)
    if (any(interactors.ok$significant)){
      interactors.genes <- paste(interactors.ok[interactors.ok$significant == TRUE, ]$gene, collapse=';')
    }
  } 
  return(interactors.genes)

}


## get all interactars in the dataset, regardless of significane
get.indata.interactors <- function(data, known.interactors){
  
  known.interactors.lst <- as.character(known.interactors[known.interactors$significant == TRUE, ]$gene)
  in.data = data[data$gene %in% known.interactors.lst, ]
  interactors.genes <- paste(in.data$gene, collapse=';')
  return(interactors.genes)

}



get.summary.aggregate <- function(directory){
  
  #

  
  file = list.files('~/Projects/03_MICOM/derived/', recursive = T, full.names = T, pattern = '17JAN.+SUMMARY')

  
  library(data.table)
  
  micom_summary <- lapply(unique(file), function(x) {
    print(x)
    read.csv(unique(x))[1,]
  })
  tab <- do.call(rbind, micom_summary)
  #tab <- tab[!duplicated(tab$file), ]
  tab$file <- as.character(tab$file)
  tab <- data.frame(tab)
  
  
  tab$martin <- gsub('03_MICOM','',tab$file)
  tab$martin <- as.numeric(gsub('martin','',regmatches(tab$martin,regexpr("M|martin[0-9]+",tab$martin))))
  tab <- tab[order(tab$martin), ]
  tab$data.bait.enriched[is.na(tab$data.bait.enriched)] <- FALSE
  
  tab$data.bait.found <- as.numeric(tab$data.bait.found)
  tab$data.bait.enriched <- as.numeric(tab$data.bait.enriched)
  
  tab1 <- tab
  tab2 <- read.csv('~/Projects/03_MICOM/14JAN2019_analysis_summary_b.csv')
  
  # by cell type
  for (mybait in unique(tab2$Bait)){
    
    all_entries <- tab2[tab2$Bait == mybait, ]
    rep <- tab2[tab2$Bait == mybait & tab2$Replicate.Correlation < 0.6,]
    only_recovered <- tab2[tab2$Bait == mybait & tab2$Bait.recovered == TRUE & tab2$Bait.enriched == FALSE & tab2$Replicate.Correlation < 0.6,]
    enriched <- tab2[tab2$Bait == mybait & tab2$Bait.recovered == TRUE & tab2$Bait.enriched == TRUE & tab2$Replicate.Correlation < 0.6,]
    
    
    write(paste('\n#######',mybait),stdout())
    print(paste('Rep less than', nrow(rep)))
    print(paste('only recov:', nrow(only_recovered)))
    print(paste('enriched and recov:', nrow(enriched)))
    
    
  }

  
  
  
  
  
  
  #tab$data.bait.enriched <- as.character(tab$data.bait.enriched)
  #tab$data.bait.enriched[tab$data.bait.enriched == 'FALSE'] <- 'BLACK'
  #tab$data.bait.enriched[tab$data.bait.enriched == 'TRUE'] <- 'RED'
  #plot(tab$rep.cor~tab$martin, col = tab$data.bait.enriched)
  
  sum(tab$data.bait.enriched)/41
  
  
}













x <- tab
x <- x[order(x$bait), ]
x <- x[order(x$cell), ]
#write.csv(x, '~/Projects/03_MICOM/analysis.summary_ordered.csv')

dat <- x[, c('bait', 'cell', 'rep.cor', 'martin', 'data.bait.enriched', 'data.bait.found')]
dat$status <- ifelse(dat$data.bait.enriched, 'Bait Significantly Enriched', 'Bait Not Found')
dat$status[(!dat$data.bait.enriched) & (dat$data.bait.found)] <- 'Bait Not Significantly Enriched'

#boxplot(dat$rep.cor~dat$bait, las=2, xlab = '', ylab = expression('R^2'), col = dat$col)

## rename to standaridize
dat$bait <- as.character(dat$bait)
dat$cell <- as.character(dat$cell)
dat$bait[dat$bait %in% c("KIAA1462(JCAD)", "JCAD (alias:KIAA1462)", "KIAA1462")] <- "JCAD (KIAA1462)"
dat$cell[dat$cell == 'teloHAEC'] <- 'Human Aortic Endothelial Cells (HAEC)'
dat$cell[dat$cell == 'EC'] <- 'Human Aortic Endothelial Cells (HAEC)'
dat$cell[dat$cell == 'SMC'] <- 'Human Coronary Artery Smooth Muscle Cells (HCASMC)'
dat$cell[dat$cell == 'THP'] <- 'THP1'
dat <- dat[order(dat$bait), ]

# Change box plot colors by groups
library(ggplot2)
plt <- ggplot(dat, aes(x=bait, y=rep.cor, color = status, shape = cell)) +
  geom_point(size = 5) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab('Replicate Correlation') + labs(title='Replicate correlation by bait and cell type.') +
  scale_color_manual(values=c("grey",'orange3',"orangered3"))
plt

ggsave('~/Projects/03_MICOM/Presentations/13Jan2020_replicate_correlation_by_bait.png', plt)
