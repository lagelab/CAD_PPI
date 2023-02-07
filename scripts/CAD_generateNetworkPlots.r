##########################################################################################
## (1) Generate network plots for PPIs identified in EC, SMC, or ALL IP-MS datasets
## (2) Generate pie charts for each network to show:
##		% of genes interacting with multiple baits
##		% of interactions in InWeb database
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

# Code tested with R version 4.2.2
library(igraph) # v1.3.5
library(qgraph) # v1.9.3
library(ggplot2) # v3.4.0
library(RColorBrewer) # v1.1.3
library(genoppi) # v1.0.13

# ----------------------------------------------------------------------------------------
# read in PPI data
allTable <- NULL
prefixes <- c('ARHGEF26_EC','BCAS3_EC','EDN1_EC','FLT1_EC','FN1_EC',
	'HDAC9_EC','JCAD_EC','PHACTR1_EC','PLPP3_EC',
	'ADAMTS7_SMC','ARHGEF26_SMC','BCAS3_SMC','EDN1_SMC','EDNRA_SMC',
	'FN1_SMC','HDAC9_SMC','JCAD_SMC','PHACTR1_SMC','PLPP3_SMC')

allTable <- read.table('../data/CAD_MasterInteractorTable.txt',
	header=T,sep='\t',stringsAsFactors=F)
allTable <- subset(allTable,
	ListName %in% prefixes & Interactor)[,c('CellType','Bait','GeneName')]


# get InWeb interactors of baits
uniqBaits <- unique(allTable$Bait)
inwebDf <- NULL
for (bait in uniqBaits) {
	tempDf <- genoppi::get_inweb_list(bait)
	if (!is.null(tempDf)) { # if bait in InWeb
		inwebInts <- subset(tempDf,significant==T)$gene
		inwebDf <- rbind(inwebDf,data.frame(Bait=bait,gene=inwebInts))
	}
}


# define rows to include in each network: EC, SMC, ALL
ecInds <- allTable$CellType=="EC"
smcInds <- allTable$CellType=="SMC"
allInds <- rep(TRUE,nrow(allTable)) 

indsList <- list(EC=ecInds,SMC=smcInds,ALL=allInds)


# iterate through each network to generate plots
for (netName in names(indsList)) {
	
	# create network
	intTable <- unique(subset(allTable,indsList[[netName]])[,c("Bait","GeneName")])
	intNetwork <- graph_from_data_frame(intTable, directed=F)

	# inInWeb boolean vector (for each edge)
	intTable$inInWeb <- NA
	for (i in 1:nrow(intTable)) {
		intTable$inInWeb[i] <- (intTable$GeneName[i] %in% 
			inwebDf$gene[inwebDf$Bait==intTable$Bait[i]])
	}

	# order columns/rows in intTable to plot bait nodes and InWeb edges last/on top
	intTable <- intTable[order(intTable$inInWeb),]
	
	# create network
	intNetwork <- graph_from_data_frame(intTable[,c('GeneName','Bait')], directed=F)

	# set vertex attributes based on # of linked baits (for node color, size, label)
	V(intNetwork)$isBait <- V(intNetwork)$name %in% unique(intTable$Bait)

	maxNumBaits <- max(table(intTable$GeneName))
	nodeColors <- brewer.pal(maxNumBaits+1,"Purples")[2:(maxNumBaits+1)]

	V(intNetwork)$size <- NA
	V(intNetwork)$color <- NA
	for (x in V(intNetwork)$name) {
		if (x %in% unique(intTable$Bait)) {
			V(intNetwork)$size[V(intNetwork)$name==x] <- maxNumBaits
			V(intNetwork)$color[V(intNetwork)$name==x] <- "red"
		} else {
			V(intNetwork)$size[V(intNetwork)$name==x] <-
				sum(intTable$GeneName==x)
			V(intNetwork)$color[V(intNetwork)$name==x] <-
				nodeColors[sum(intTable$GeneName==x)]
		}
	}

	V(intNetwork)$label <- NA
	V(intNetwork)$label[V(intNetwork)$size>=2] <-
		V(intNetwork)$name[V(intNetwork)$size>=2]

	# set edge attributes based on InWeb data (for edge color)
	E(intNetwork)$inInWeb <- intTable$inInWeb

	# network plot layout
	e <- get.edgelist(intNetwork,names=F)
	l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(intNetwork),niter=5000,
		area=(vcount(intNetwork)^1.7),repulse.rad=(vcount(intNetwork)^2.2))
								
	# generate network plot
	plotName <- paste("../output/CAD_",netName,"_NetworkPlot.pdf",sep="")
	pdf(plotName,width=7,height=7)

	# with gene name labels (for baits + interactors linked to > 1 bait)
	print(plot(intNetwork, edge.width=0.5,
		vertex.color=V(intNetwork)$color,
		vertex.size=V(intNetwork)$size+1, vertex.frame.color="grey30",
		vertex.label=V(intNetwork)$label, vertex.label.color="black",
		vertex.label.cex=0.08,
		edge.color=c("grey","blue")[1+(E(intNetwork)$inInWeb)],
		layout=l))

	# no gene name labels
	print(plot(intNetwork, edge.width=0.5,
		vertex.color=V(intNetwork)$color,
		vertex.size=V(intNetwork)$size+1, vertex.frame.color="grey30",
		vertex.label=NA, vertex.label.cex=0.08,
		edge.color=c("grey","blue")[1+(E(intNetwork)$inInWeb)],
		layout=l))

	dev.off()


	# --------------------------------------------------------------------------------
	# generate pie charts
	pieName <- paste("../output/CAD_",netName,"_PieCharts.pdf",sep="")
	pdf(pieName,width=3,height=3)

	# % of genes that are interactors of diff. # of baits
	genePieTable <- data.frame(table(table(intTable$GeneName)))
	colnames(genePieTable) <- c("NumBaits","GeneCount")
	genePieTable$Percent <- genePieTable$GeneCount/sum(genePieTable$GeneCount)*100

	print(ggplot(genePieTable, aes(x="", y=GeneCount, fill=NumBaits)) +
	geom_bar(width=1, stat="identity", color="black") +
	coord_polar("y") + 
	scale_fill_manual(name="# of Baits", values=nodeColors) +
	theme_void())

	# % of interactions in vs. not in InWeb
	inwebTable <- data.frame(table(intTable$inInWeb))
	inwebTable$Percent <- inwebTable$Freq/sum(inwebTable$Freq)*100
	names(inwebTable)[1] <- 'IsInWeb'
		
	print(ggplot(inwebTable, aes(x="", y=Freq, fill=IsInWeb)) +
	geom_bar(width=1, stat="identity", color="black") +
	coord_polar("y") + scale_fill_manual(name="InWeb", values=c("grey","blue")) +
	theme_void())

	dev.off()

	print(genePieTable)

	print(inwebTable)


	# --------------------------------------------------------------------------------
	# gene frequency table (for identifying top recurrent interactors)
	geneFreqTable <- NULL
	for (i in 1:nrow(intTable)) {
		isGene <- intTable$GeneName==intTable$GeneName[i]
		geneFreqTable <- rbind(geneFreqTable,
			data.frame(Interactor=intTable$GeneName[i],NumBaits=sum(isGene),
			Baits=paste(intTable$Bait[isGene],collapse=", ")))
	}	
	geneFreqTable <- unique(geneFreqTable)

	# order by decreasing frequency, then by alphabetical
	geneFreqTable <- geneFreqTable[order(-geneFreqTable$NumBaits,
		as.character(geneFreqTable$Interactor)),]

	# output to file
	freqName <- paste("../output/CAD_",netName,"_InteractorFrequency.txt",sep="")
	write.table(geneFreqTable,freqName,sep="\t",row.names=F,quote=F)
}

