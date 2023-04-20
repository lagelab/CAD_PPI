


ggbarplot <- function(df, bonf, colors = c("red","sienna1","black"), ylab = 'Tissue type'){
  
  stopifnot(c('list_name','pvalue') %in% colnames(df))
  
  # log-pvalue
  df$logpvalue <- -log10(df$pvalue)
  
  # set significance thresholds
  df$significance <- 'Not significant'
  df$significance[df$pvalue < 0.05] <- 'Nominally significant'
  df$significance[df$pvalue < bonf] <- 'Bonferroni threshold'
  df$significance <- as.factor(df$significance)
  
  # generate color palette
  names(colors) <- c('Bonferroni threshold','Nominally significant', 'Not significant')
  color_scale <- scale_fill_manual(name = "significance", values = colors)
  
  # generate plot
  plt <- ggplot(df, aes(y = reorder(list_name, logpvalue), x = logpvalue, fill = significance)) + 
    geom_bar(stat = 'identity', alpha = 0.9) +
    geom_vline(xintercept = -log10(0.05), color = "black", linetype = 'dashed') + 
    geom_vline(xintercept = -log10(bonf), color = "black", linetype = 'dashed') + 
    xlab('-log10(P-value)') + ylab(ylab) +
    color_scale +
    theme_bw() 
  
  return(plt)
  
}
