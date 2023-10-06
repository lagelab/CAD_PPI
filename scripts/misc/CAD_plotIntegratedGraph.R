
## Make graph
# generate network plots
# generate pie charts:
#	(1) % of genes interacting with multiple baits
#	(2) % of interactions in InWeb

rm(list=ls())
#setwd("~/Google Drive/LageLab/T2D")

library(igraph)
library(qgraph)
library(hash)
library(ggplot2)
library(RColorBrewer)

# InWeb interactors
#load("~/Google Drive/LageLab/BINe_IPs/Genoppi-master/data/InWeb_combined_Oct2018.RData")

# what baits do we have?
info = read.csv('10FEB2020_masterlist_analysis_plan.csv')
baits = c(as.character(unique(info$Bait)), 'ARHGEF26-MT')

# gather integrated plots
#files <- list.files("data/interactor_lists/run056/", full.names = T)
#files <- files[!grepl(pattern = '^(COMBINED)',files)]
#lapply(files, function(x){
  fread)x
})


# interactor table
all_table <- list(
  #read.table("data/interactor_lists/run056/COMBINED_EC-SMC.InteractorTable.txt",header=T,sep="\t",stringsAsFactors=F),
  read.table("data/interactor_lists/run056/COMBINED_SMC.InteractorTable.txt",header=T,sep="\t",stringsAsFactors=F),
  read.table("data/interactor_lists/run056/COMBINED_EC.InteractorTable.txt",header=T,sep="\t",stringsAsFactors=F)
)

all_table <- do.call(rbind, all_table)
all_table$significant <- all_table$Interactor
all_table$cell <- all_table$CellType
all_table$gene <- all_table$GeneName
all_table$bait <- all_table$Bait
all_table$Bait

#all_table <- read.table('06MAR2020_micom_tier_1_ips.csv')

## genetics
gene_table = read.table('derived/genelists/cad_gwas_genes_prio1_2.txt')
colnames(gene_table) = c('gene')

# network names
networkList = c('SMC','EC')
#networkList <- c("Tier1_BetaCell","Tier1_Hepatocyte","Tier1_Combined",
#                 "Tier2_BetaCell","Tier2_Hepatocyte","Tier2_Adipocyte","Tier2_Combined")
#tierList <- c(1,1,1,2,2,2,2)
cell_list <- c('SMC','EC')
graphics.off()

for (cell in cell_list) {

  # filter full interactor table for interactors in this specfic network
  #isCell <- rep(TRUE,nrow(allTable))
  #
  #if (!is.na(cellList[n])) { isCell <- allTable$CellType==cellList[n]}
  #df_int <- allTable[allTable$Tier<=tierList[n] & isCell & allTable$Interactor=="y",
  #                     c("Bait","GeneName")]
  df_int = all_table[all_table$cell %in% cell & all_table$significant == T, ]
  df_int = df_int[df_int$gene %nin% baits, ]
  df_int <- unique(df_int)

  # create network (+ inInWeb boolean vector)
  #inInWeb <- c()
  #networkEdges <- c()
  #for (i in 1:nrow(df_int)) {
  #  inInWeb <- c(inInWeb, df_int$gene[i] %in% inweb_hash[[df_int$Bait[i]]])
  #  networkEdges <- c(networkEdges, df_int$Bait[i], df_int$gene[i])
  #}
  
  tmp = unique(na.omit(c(as.character(df_int$bait), as.character(df_int$gene))))
  nodes = data.frame(name=tmp, id = 1:length(tmp))
  edges = data.frame(source=df_int$bait, target = df_int$gene, weight = 1)
  edges = edges[complete.cases(edges),]
  nodes$name <- as.character(nodes$name)
  edges$source <- as.character(edges$source)
  edges$target <- as.character(edges$target)
  int_network <- graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)
  
  #in_inweb = df_int$known_inweb_interactor
  #network_edges = df_int[,c('bait','gene')]
  #int_network <- graph(network_edges,directed=F)

  # set vertex attributes (for node color, size, label)
  V(int_network)$bait <- V(int_network)$name %in% unique(df_int$bait)
  max_num_baits <- max(table(df_int$gene))
  node_colors <- brewer.pal(max_num_baits+2,"Purples")[2:(max_num_baits+2)]

  V(int_network)$size <- NA
  V(int_network)$color <- NA
  for (x in V(int_network)$name) {
    if (x %in% unique(df_int$bait)) {
      V(int_network)$size[V(int_network)$name==x] <- max_num_baits
      V(int_network)$color[V(int_network)$name==x] <- "red"
    } else {
      V(int_network)$size[V(int_network)$name==x] <- sum(df_int$gene==x)
      V(int_network)$color[V(int_network)$name==x] <- node_colors[sum(df_int$gene==x)]
    }
    if (x %in% gene_table$gene & x %nin% baits){
      V(int_network)$color[V(int_network)$name==x] <- "orange"
    }
  }

  V(int_network)$label <- NA
  V(int_network)$label[V(int_network)$size>=2] <- V(int_network)$name[V(int_network)$size>=2]

  # set edge attributes (for edge color)
  E(int_network)$inInWeb <- df_int$known_inweb_interactor
  
  # calculate hypergeometric overlap
  baits_found = sum(baits %in% gene_table$gene)
  N = length(unique(all_table$gene))  # Size of population
  n = length(gene_table$gene)-baits_found # number of success in population
  a = length(V(int_network)$name) # draws
  x = sum(V(int_network)$name %in% gene_table$gene)-baits_found # found overlap
  pvalue = dhyper(sum(V(int_network)$name %in% gene_table$gene),m=a,n=N-a,k=n)

  
  # network plot layout
  e <- get.edgelist(int_network,names=F)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(int_network),niter=5000,
                                         area=(vcount(int_network)^1.7),repulse.rad=(vcount(int_network)^2.2))
  
  
  # ----------------------------------------------------------------------------------------
  # generate network plot
  pdf(paste("derived/plots/210406",cell,"NetworkPlot.pdf",sep="_"),width=7,height=7)

  plot(int_network, edge.width=0.5,
       vertex.color=V(int_network)$color,
       vertex.size=V(int_network)$size+1, vertex.frame.color="grey30",
       vertex.label=V(int_network)$label, vertex.label.color="black",
       vertex.label.cex=0.08,
       edge.color=c("grey","blue")[1+(E(int_network)$inInWeb)],
       layout=l)
  legend("bottomright", paste('CAD (Genetics) P-value:',pvalue), bty="n", cex = 0.7)
  
  plot(int_network, edge.width=0.5,
       vertex.color=V(int_network)$color,
       vertex.size=V(int_network)$size+1, vertex.frame.color="grey30",
       vertex.label=NA, vertex.label.cex=0.08,
       edge.color=c("grey","blue")[1+(E(int_network)$inInWeb)],
       layout=l)
  
  plot(0:100,dhyper(x=0:100,m=a,n=N-a,k=n),type = 'h',main = 'probability of x succeses')
  abline(v=x, col = 'red'); 


                
  # ----------------------------------------------------------------------------------------
  # generate table of node connections
  
  int_network_names = names(sort(degree(int_network), decreasing = T))
  int_network_degree = sort(degree(int_network), decreasing = T )
  int_dat = data.frame(gene = int_network_names, degree = int_network_degree)
  rownames(int_dat) = NULL
  #write.table(int_dat,paste("tables/06MAR2020",cell,"node_connections.txt",sep="_"), quote = F, sep = '\t', row.names = F)

  # ----------------------------------------------------------------------------------------
  # generate pie charts

  pdf(paste("plots/210406",cell,"PieCharts.pdf",sep="_"),width=3,height=3)

  # % of genes that are interactors of diff. # of baits
  genePieTable <- data.frame(table(table(df_int$gene)))
  colnames(genePieTable) <- c("NumBaits","GeneCount")
  genePieTable$Percent <- genePieTable$GeneCount/sum(genePieTable$GeneCount)*100

  #genePieTable$NumBaits <- factor(genePieTable$NumBaits,levels=rev(genePieTable$NumBaits))
  #y.breaks <- cumsum(rev(genePieTable$GeneCount)) - rev(genePieTable$GeneCount)/2
 
  print(ggplot(genePieTable, aes(x="", y=GeneCount, fill=NumBaits)) +
          geom_bar(width=1, stat="identity", color="black") +
          coord_polar("y") + scale_fill_manual(name="# of Baits", values=node_colors) + #palette="Purples") +
          theme_void())
  #+ theme(axis.text.x=element_text(color="black"))
  #scale_y_continuous(breaks=y.breaks, labels=rev(genePieTable$Percent))

  # % of interactions in vs. not in InWeb
  inwebTable <- data.frame(table(df_int$known_inweb_interactor))
  inwebTable$Percent <- inwebTable$Freq/sum(inwebTable$Freq)*100
  
  print(ggplot(inwebTable, aes(x="", y=Freq, fill=Var1)) +
          geom_bar(width=1, stat="identity", color="black") +
          coord_polar("y") + scale_fill_manual(name="InWeb", values=c("grey","blue")) +
          theme_void())

  dev.off()
  
  print(cell)
  print(genePieTable)
  print(inwebTable)

}
graphics.off()
