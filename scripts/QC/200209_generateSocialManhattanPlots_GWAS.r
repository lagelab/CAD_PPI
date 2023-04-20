
rm(list=ls())

#.libPaths('~/Toolbox/rlib')
setwd('~/Projects/03_MICOM/')
library(ggplot2)
library(dplyr)
library(ggrepel)


# read in genetic data

## load anneals
anneals <- read.table('derived/genelists/27AUG20_haarst2017_ensembl_anneal_mapping.tsv', sep = '\t', stringsAsFactors = F, header = T)

anneals <- anneals[anneals$gene_type == 'protein_coding',]
mergeDf <- anneals
mergeDf$P <- mergeDf$P.value

## check presence of certain genes
y1 = unlist(strsplit(y, split = '(\\ )|(\n)'))
y2 = y1[y1 !='']

# rename columns to match yu-hans format
colnames(mergeDf)[colnames(mergeDf) == 'chromosome_name'] <- 'Chr'
colnames(mergeDf)[colnames(mergeDf) == 'geneStart'] <- 'Pos'
colnames(mergeDf)[colnames(mergeDf) == 'external_gene_name'] <- 'Gene'

# remove columns 
mergeDf <- mergeDf[,c('Gene', 'SNP', 'P', 'Chr', 'Pos')]

# deal with synonyms
c('KIAA1462', 'JCAD') %in% mergeDf$Gene
c('PLPP3', 'PPAP2B') %in% mergeDf$Gene
mergeDf$Gene[mergeDf$Gene == 'KIAA1462'] <- 'JCAD'
mergeDf$Gene[mergeDf$Gene == 'PPAP2B'] <- 'PLPP3'

# set minimum P to 1e-30
mergeDf$P[mergeDf$P<1e-30] <- 1e-30
mergeDf <- mergeDf[order(mergeDf$Chr,mergeDf$Pos),]
dim(mergeDf)



nChr <- length(unique(mergeDf$Chr))
mergeDf$Coord <- NA
s <- 0
#nbp <- c()
nbp <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,
	141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,
	81195210,78077248,59128983,63025520,48129895,51304566) # GRCh37 chr 1-22 lengths
chrTicks <- NULL
for (i in 1:22) { #unique(mergeDf$Chr)) {
	#nbp[i] <- max(mergeDf[mergeDf$Chr==i,"Pos"])
	mergeDf[mergeDf$Chr==i,"Coord"] <- mergeDf[mergeDf$Chr==i,"Pos"] + s
	chrTicks <- c(chrTicks,floor(nbp[i]/2) + s)
	s <- s + nbp[i]
	#chrTicks <- c(chrTicks,mergeDf[mergeDf$Chr==i,"Coord"][floor(length(mergeDf[mergeDf$Chr==i,"Coord"])/2)+1])
}



# generate and combine network lists

combine_int <- function(f){
  df = read.csv(f, sep = '\t')
  df = df[ ,c('Bait', 'GeneName', 'Interactor', 'CellType')]
  colnames(df) <- c('Bait','Gene','Significant','Cell')
  return(df)
}


## table for EC 
path = 'derived/maps/maps-run056/'
filesEC <- list.files(path,'EC\\.InteractorTable\\.txt', full.names = T)
filesEC <- filesEC[!grepl('COMBINED',filesEC)]
tabEC <- do.call(rbind, lapply(filesEC, combine_int))

## Table for SMC
path = 'derived/maps/maps-run056/'
filesSMC <- list.files(path,'SMC\\.InteractorTable\\.txt', full.names = T)
filesSMC <- filesSMC[!grepl('COMBINED',filesSMC)]
tabSMC<- do.call(rbind, lapply(filesSMC, combine_int))

# combine network lists
#tabEC <- read.csv('derived/maps/11FEB2020_EC_map.csv',header=T,stringsAsFactors=F)[,c(2,3,4)]
#tabHEK <- read.csv('derived/maps/11FEB2020_HEK293_map.csv',header=T,stringsAsFactors=F)[,c(2,3,4)]
#tabSMC <- read.csv('derived/maps/11FEB2020_SMC_map.csv',header=T,stringsAsFactors=F)[,c(2,3,4)]
#tabEC$cell <- 'EC'
#tabSMC$cell <- 'SMC'
#tabHEK$cell <- 'HEK'
allTable <- rbind(tabEC,tabSMC)
allTable$Tier <- 1
colnames(allTable) <- c('Bait', 'Gene', 'Significant', 'Cell', 'Tier')

## only keep significant data
nrow(allTable)
allTable <- allTable[allTable$Significant == TRUE, ]
nrow(allTable)

# check what baits are in network
mybaits = unique(allTable$Bait)
mybaits[!mybaits %in% mergeDf$Gene]

### deal with EDN1 (experimentaly linked to PHACTR1)
if (T){
  rowEDN1 <- mergeDf[mergeDf$Gene == 'PHACTR1',] 
  rowEDN1$Gene <- 'EDN1'
  rowEDN1$Pos <- 12290596
  mergeDf <- rbind(rowEDN1, mergeDf)
  ### 
}


# exclude "FN1" "PLPP3" "ARHGEF26(WT)" "BCAS3" b/c bait not in gwas data (but are in exome chip data)
allTable <- allTable[allTable$Bait %in% mergeDf$Gene,]  
networkList <- c('SMC_analysis', 'EC_analysis')
tierList <- c(1, 1)
cellList <- c('SMC', 'EC')

# iterate through each network
for (n in 1:length(networkList)) {

	# filter full interactor table for interactors in this specfic network
	isCell <- rep(TRUE,nrow(allTable))
	if (!is.na(cellList[n])) { isCell <- allTable$Cell==cellList[n]}
	intTable <- allTable[allTable$Tier<=tierList[n] & isCell & allTable$Significant==T,c("Bait","Gene")]
	intTable <- unique(intTable)

	# create data frame for genes that have P < 1e-5 and are baits or interactors in network
	plotDf <- merge(mergeDf,intTable,by='Gene')
	plotDf$Bait_P <- NA
	plotDf$Bait_Coord <- NA
	
	baitList <- unique(intTable$Bait)
	for (Bait in baitList) {
		# add bait coordinates (for plotting social links between genes and baits)
	  plotDf[plotDf$Bait==Bait,"Bait_P"] <- mergeDf$P[mergeDf$Gene==Bait]
		plotDf[plotDf$Bait==Bait,"Bait_Coord"] <- mergeDf$Coord[mergeDf$Gene==Bait]
	}

	for (Bait in baitList) {
	  print(Bait)
		baitDf <- data.frame(mergeDf[mergeDf$Gene==Bait,c("Gene","SNP","P","Chr","Pos","Coord")],
		Bait=NA,Bait_P=NA,Bait_Coord=NA)
		plotDf <- rbind(plotDf,baitDf,row.names=NULL)
	}
	
	# write out file
	filename = paste("plots/210303",networkList[n],"SocialManhattanPlotData_GWAS_EDN1.csv",sep="_")
	write.csv(plotDf[order(plotDf$Bait),], filename, quote = F, row.names = F)

	dim(plotDf)
	plotDf$text <- paste0(apply(plotDf[,c('Gene','SNP')], 1, paste, collapse = ' ('), ')')

	# significance group (for node color)
	plotDf[plotDf$P<1e-5 | plotDf$Gene %in% baitList,]	
	plotDf$Significance <- 'FDR < 5%'
	#plotDf[plotDf$P < e-4,"Significance"] <- "P < 1e-4"
	#plotDf[plotDf$P < 1e-5,"Significance"] <- "P < 1e-5"
	plotDf[plotDf$P < 1e-5,"Significance"] <- "P < 1e-5" # "loci-wide significant"
	plotDf[plotDf$P < 5e-8,"Significance"] <- "P < 5e-8"
	plotDf[plotDf$Gene %in% baitList,"Significance"] <- "BAIT"
	plotDf$Significance <- factor(plotDf$Significance,levels=c("BAIT","P < 5e-8","P < 1e-5", 'FDR < 5%'))

	# other plotting parameters
	ylim <- abs(floor(log10(min(plotDf$P)))) + 0.5 
	sig <- 5e-8 

	uniqDf <- unique(plotDf[,c("Gene","Coord","P","Significance",'text')])

	# for now remove, annotation og sub-GWS
	uniqDf$Significance <- as.character(uniqDf$Significance)
	uniqDf$Significance[uniqDf$Significance != 'BAIT'] <- 'INTERACTOR'
	uniqDf$Significance <- as.factor(uniqDf$Significance)
	
	subt = 'Genelist derived from 161 (158) CAD SNP anneals (r2 > 0.6 and 50kb upstream/downstream)'
	
	# Manhattan plot
	pdf(paste("plots/210303",networkList[n],"SocialManhattanPlot_GWAS_EDN1.pdf",sep="_"), width=11,height=7)

	print(n)
	print(plotDf)
	print(uniqDf)
	
	# P < 1e-5 version
	print(ggplot(uniqDf,mapping=aes(x=Coord,y=-log10(P),color=as.factor(Significance))) +
		#geom_hline(yintercept=-log10(sig),color="orange",linetype="dashed") +
		
		geom_segment(plotDf,mapping=aes(x=Coord,y=-log10(P),
			xend=Bait_Coord,yend=-log10(Bait_P)),color="grey") +

		geom_point(size=2) +
		geom_text_repel(mapping=aes(label=Gene),size=2.5,show.legend=FALSE) +

		scale_color_manual(name="MICOM Gene Groups",#labels=c("Baits for IP","P < 5e-8","P < 1e-5"),
			values=c("red","blue","darkcyan", 'cyan')) +
			
		scale_x_continuous(limits=c(0,max(mergeDf$Coord)),label=c(as.character(1:22)),breaks=chrTicks) +
		scale_y_continuous(expand=c(0,0), limits=c(0, ylim)) +

		ggtitle("Social Manhattan Plot (GWAS, Pim van der Haarst 2017)", subt) +
		xlab("Chromosomal Position") + ylab("-log10(P-value)") +
		theme_classic())

	# P < 1e-5 version
	print(ggplot(uniqDf,mapping=aes(x=Coord,y=-log10(P),color=as.factor(Significance))) +
	        #geom_hline(yintercept=-log10(sig),color="orange",linetype="dashed") +
	        
	        geom_segment(plotDf,mapping=aes(x=Coord,y=-log10(P),
	                                        xend=Bait_Coord,yend=-log10(Bait_P)),color="grey") +
	        
	        geom_point(size=2) +
	        geom_text_repel(mapping=aes(label=text),size=2.5,show.legend=FALSE) +
	        
	        scale_color_manual(name="MICOM Gene Groups",#labels=c("Baits for IP","P < 5e-8","P < 1e-5"),
	                           values=c("red","blue","darkcyan", 'cyan')) +
	        
	        scale_x_continuous(limits=c(0,max(mergeDf$Coord)),label=c(as.character(1:22)),breaks=chrTicks) +
	        scale_y_continuous(expand=c(0,0), limits=c(0, ylim)) +
	        
	        ggtitle("Social Manhattan Plot (GWAS, Pim van der Haarst 2017)", subt) +
	        xlab("Chromosomal Position") + ylab("-log10(P-value)") +
	        theme_classic())
	
	
	# P < 5e-8 version
	#print(ggplot(uniqDf[uniqDf$P<5e-8,],mapping=aes(x=Coord,y=-log10(P),color=as.factor(Significance))) +
	#	#geom_hline(yintercept=-log10(sig),color="orange",linetype="dashed") +
#		
	#	geom_segment(plotDf[plotDf$P<5e-8,],mapping=aes(x=Coord,y=-log10(P),
#			xend=Bait_Coord,yend=-log10(Bait_P)),color="grey") +
#
#		geom_point(size=2) +
#		geom_text_repel(mapping=aes(label=Gene),size=2.5,show.legend=FALSE) +
#
#		scale_color_manual(name="MICOM Gene Groups",#labels=c("Baits for IP","P < 5e-8","P < 1e-5"),
#			values=c("red","blue","darkcyan", 'cyan')) +
#		
#		scale_x_continuous(limits=c(0,max(mergeDf$Coord)),label=c(as.character(1:22)),breaks=chrTicks) +
#		scale_y_continuous(expand=c(0,0), limits=c(0, ylim)) +
#		
#		ggtitle("Social Manhattan Plot (GWAS, Pim van der Haarst 2017)", subt) +
#		xlab("Chromosomal Position") + ylab("-log10(P-value)") +
#		theme_classic())

	# baits only version
	print(ggplot(uniqDf[uniqDf$Gene %in% baitList,],mapping=aes(x=Coord,y=-log10(P),color=as.factor(Significance))) +
		#geom_hline(yintercept=-log10(sig),color="orange",linetype="dashed") +
		
		geom_point(size=2) +
		geom_text_repel(mapping=aes(label=Gene),size=2.5,show.legend=FALSE) +

		scale_color_manual(name="MICOM Gene Groups",#labels=c("Baits for IP","P < 5e-8","P < 1e-5"),
			values=c("red","blue","darkcyan", 'cyan')) +
		
		scale_x_continuous(limits=c(0,max(mergeDf$Coord)),label=c(as.character(1:22)),breaks=chrTicks) +
		scale_y_continuous(expand=c(0,0), limits=c(0, ylim)) +
		
		ggtitle("Social Manhattan Plot (GWAS, Pim van der Haarst 2017)") +
		xlab("Chromosomal Position") + ylab("-log10(P-value)") +
		theme_classic())

	dev.off()

}

