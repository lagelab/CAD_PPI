setwd('~/Projects/03_MICOM/')

## deprecated (not used anymore)


library(ggplot2)

# ready paths and prefixes
date = '201120'
run = 'run056'
target = 'derived/maps/maps-run056/'

# setup paths and run data
paths = read.csv('run056_filtered_tier1_paths.tsv', sep = '\t')$x[1]

pdf('~/Desktop/201129_inweb_test.pdf')
for (p in paths){
  
  print(p)
  #p = "data/genoppi_input/run056/Whitehead.m10.JCADvsMock.SMC.3xFLAG.20NOV2020.tsv" 
  df <- read.csv(p, sep = '\t')
  df$gene[df$gene == 'JCAD'] <- 'KIAA1462'
  
  # save inweb edges to memory
  inweb <- lapply(unique(df$gene), function(x) {d <- get_inweb_list(x); print(x); d$gene[d$significant]})
  names(inweb) <- unique(df$gene)
  
  # count edges
  edges <- function(query, data, inweb){
    
    if (query %in% data$gene){
      if (query %in% names(inweb)){
        
        count = sum(inweb[[query]] %in% df$gene)
        known = length(inweb[[query]])
        sig = sum(inweb[[query]] %in% df$gene[df$significant])
        notsig = sum(inweb[[query]] %in% df$gene[!df$significant])
        return(data.frame(n = count, known = known, n.sig = sig, n.notsig = notsig))
        
      }
    }
  }
  
  
  #res <- do.call(rbind, lapply(df$gene, function(x){
  #  return(edges(x, df, inweb))
  #}))
  #df <- cbind(df, res)
  
  #p1 = ggplot(df, aes(y=-log10(pvalue), x = logFC, color = significant, size = n/known)) +
  #  ylab('Inweb Interactors / Known InWeb interactors') +
  #  geom_point() +
  #  ggtitle(basename(p)) +
  #  theme_minimal()
  
  #p2 = ggplot(df, aes(y=n/known, x = -log10(FDR), color = significant)) +
  #  ylab('Inweb Interactors / Known InWeb interactors') +
  #  geom_point() +
  #  ggtitle(basename(p)) +
  #  theme_minimal()
  
  #p3 = ggplot(df, aes(y=n/known, x = significant, fill = significant)) +
  #  ylab('Inweb Interactors / Known InWeb interactors') +
  #  geom_boxplot() +
  #  geom_jitter() +
  #  ggtitle(basename(p)) +
  #  theme_minimal()
  
  #print(p1)
  #print(p2)
  #print(p3)
  
  
  ## permute data
  #set.seed(4231)
  #dfp <- df
  #indicies <- sample(1:nrow(df), nrow(df))
  #dfp$gene <- dfp$gene[indicies]
  
  #dfp$edge <- unlist(lapply(dfp$gene, function(x){
  #  return(edges(x, df, inweb))
  #}))
  
  #ggplot(dfp, aes(y=edge, x = significant)) +
  #  geom_boxplot()
  
}
graphics.off()
