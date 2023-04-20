# Get a data.frame witn the filename, replicate correlation
# whether bait was recovered, and whether the file contains
# mock
get_micom_ip_summary <- function(data){
  
  # find successfull IPs
  bait_recovered = unlist(lapply(names(data), function(dname){
    header = unlist(strsplit(dname,'\\.'))
    bait = gsub('\\((W|M)T\\)','',unlist(strsplit(header[2], 'vs'))[1], perl = T)
    df = data[[dname]]
    bait_found = any(grepl(toupper(bait), toupper(df$gene)) | grepl('FLAG', toupper(df$gene)))
  }))
  
  # select ones with higher correlations
  good_rep_ip = unlist(lapply(names(data), function(dname){
    df = data[[dname]]
    corr = cor(df$rep1, df$rep2)
    return(corr)
  }))
  
  # select good IPs and remove non-IgG comparisons
  #good_data = data[good_rep_ip & bait_recovered & grepl('Mock',names(data))]
  return(data.frame(replicate_correlation=good_rep_ip, bait_found=bait_recovered, mock=grepl('Mock',names(data)), fname = names(data)))
}

# Take all genoppi IPs are prepare them in a 2xn data.frame
# in which the first column is the bait, and the 2nd column
# is the interactor
melt_ips <- function(data, inweb = T){
  interactors_overview = lapply(names(data), function(dname){
    
    # Ready data for subsetting
    inweb_bait = ''
    header = unlist(strsplit(dname,'\\.'))
    bait = unlist(strsplit(header[2], 'vs'))[1]
    bait = gsub('ARHGEF26\\(WT\\)','ARHGEF26',bait)
    bait = gsub('ARHGEF26\\(MT\\)','ARHGEF26-MT',bait)
    
    # subset data and get relevant cols
    df = data[[dname]]
    ndf = df[,c('gene','significant')]
    ndf$bait = bait
    ndf = ndf[c('bait','gene','significant')]
    ndf$cell = gsub('^HEK293T$', 'HEK293', as.character(header[3]))
    ndf$known_inweb_interactor = FALSE
    
    # get inweb interactors
    if (!is.null(bait) & inweb){
      if (inweb_bait != bait){
        inweb_bait = gsub('ARHGEF26-MT','ARHGEF26',bait)
        inweb_bait = gsub('JCAD', 'KIAA1462', inweb_bait)
        inweb_bait = gsub('PLPP3', 'PPAP2B', inweb_bait)
        inweb_interactors = interactors(inweb_bait, verbose = T)
      }
      if (any(ndf$gene %in% inweb_interactors[inweb_interactors$significant,]$gene)){
        ndf[ndf$gene %in% inweb_interactors[inweb_interactors$significant, ]$gene, ]$known_inweb_interactor = TRUE
      }
    }
    return(ndf)
  })
  return(as.data.frame(do.call(rbind, interactors_overview)))
}




#write.csv(df_int,'derived/maps/11FEB2020_non_cell_specific_map.csv')

#pdf('12FEB2020_integrated_plots_uniprot.pdf', width = 15, height = 15)
# select cell-specific maps
#for (cell in c("EC", "SMC", "HEK293")){
#  good_data_cell = data[good_rep_ip & bait_recovered & grepl('Mock',names(data)) & grepl(cell,names(data))]
#  interactors_cell_overview = lapply(names(good_data_cell), function(dname){
#    header = unlist(strsplit(dname,'\\.'))
#    bait = unlist(strsplit(header[2], 'vs'))[1]
#    #bait = gsub('\\((W|M)T\\)','',bait)
#    df = good_data_cell[[dname]]
#    ndf = df[,c('gene','significant')]
#    ndf$bait = bait
#    ndf = ndf[c('bait','gene','significant')]
#  })
#  df_int = as.data.frame(do.call(rbind, interactors_cell_overview))
#  #write.csv(df_int,paste0('derived/maps/11FEB2020_',cell,'_map.csv'))
#  g = graph_micom(df_int, cell, 12)
#  head(sort(degree(g), decreasing = T), n = 25)
#  Sys.sleep(5)
#}
#graphics.off()

# plot maps


# simple way of setting up a graph
graph_micom <- function(df_int, note, seed = 16){
  
  set.seed(seed)
  df_int = df_int[df_int$significant == TRUE,]
  tmp = unique(na.omit(c(as.character(df_int$bait), as.character(df_int$gene))))
  nodes = data.frame(name=tmp, id = 1:length(tmp))
  edges = data.frame(source=df_int$bait[df_int$significant], target = df_int$gene[df_int$significant], weight = 1)
  edges = edges[complete.cases(edges),]
  nodes$name <- as.character(nodes$name)
  edges$source <- as.character(edges$source)
  edges$target <- as.character(edges$target)
  g <- graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)
  
  # coloring
  baits = c(as.character(info$Bait), 'ARHGEF26(MT)', 'ARHGEF26(WT)')
  V(g)$color <- ifelse(V(g)$name %in% baits, "orange", "grey")
  V(g)$vertex_degree  <- log2(degree(g))
  l <- layout_with_fr(g, niter=1000)
  par(mar=c(1,1,1,1))
  plot(g, # change color of nodes
       vertex.label.color = "black", # change color of labels
       vertex.label.cex = .45, # change size of labels to 75% of original size
       edge.curved=.10, # add a 25% curve to the edges
       edge.color="grey20",
       vertex.size = V(g)$vertex_degree,
       layout = l) # change edge color to grey
  title(paste("MICOM cell-specific integrated plot (",note,")"),cex.main=1,col.main="black")
  return(g)
}


