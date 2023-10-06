# heatmap
make_cell_heatmap <- function(df, fdr.threshold = 0.5, ylab = 'Gene Ontology (MF)', main = 'Title is here.'){
  
  stopifnot(is.list(df))
  stopifnot(!is.data.frame(df))
  
  # convert into new unified data.frame 
  hm_data <- do.call(rbind, lapply(names(df), function(name){
    x = df[[name]]
    data.frame(go = x$dataset,  pvalue = x$pvalue, FDR = x$FDR, experiment = name,cell = unlist(strsplit(name, split = '_'))[2])
  }))
  
  # subset heatmap, so that only top most enriched things are plotted
  exclude <- unlist(lapply(unique(hm_data$go), function(x) all(hm_data$FDR[hm_data$go %in% x] > fdr.threshold)))
  hm_data <- hm_data[hm_data$go %in% unique(hm_data$go)[!exclude],]
  hm_data$label <- ifelse(hm_data$pvalue < 0.05,ifelse(hm_data$FDR < 0.05, '**', '*'),'')
  
  # make plots for SMC and EC cell
  combined_plots <- lapply(c("SMC",'EC'), function(x){
    p <- ggplot(hm_data[hm_data$cell %in% x,], aes(y = reorder(go, -log10(pvalue)), x = experiment, fill = -log10(pvalue), label = label)) +
      geom_tile() + theme_bw() + 
      geom_text() +
      scale_fill_gradient(low = 'white', high = ifelse(x == 'SMC', 'red','blue')) +
      xlab('IP') + 
      ylab(ylab) +
      theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
      facet_wrap(~cell) +
      theme(legend.position= unlist(ifelse(x == 'SMC', 'left','right'))) +
      if (x == "SMC") scale_y_discrete(position = "left") else scale_y_discrete(position = "right")
    return(p)
  })
  
  # prepare plot
  title <- ggdraw() + draw_label(main, fontface='bold')
  plot_row <- plot_grid(
    combined_plots[[1]],
    combined_plots[[2]]
  )
  
  # plot grid
  plot_grid(
    title, 
    plot_row,
    ncol = 1,
    rel_heights=c(0.1, 1),
    labels = NULL
  )
  
}