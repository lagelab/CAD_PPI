##########################################################################################
## Compare pairs of IP-MS datasets
## (1) Calculate pairwise similarity metrics for each pair of IPs
## (2) Generate Venn diagrams of interactors for IPs of the same bait
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

# Code tested with R version 4.2.2
library(reshape2) # v1.4.4
library(dplyr) # v1.1.0
library(ggplot2) # v3.4.0
library(ggVennDiagram) # v1.2.2

# ----------------------------------------------------------------------------------------
# read in data for 20 IP-MS datasets

summaryDf <- read.csv('../data/CAD_IpSummary.txt',header=T,sep='\t',stringsAsFactors=F)
names(summaryDf)[1] <- 'IP'

# read in data for each IP
allDf <- NULL
for (i in 1:nrow(summaryDf)) {
	
	ipDf <- read.csv(summaryDf$File.path[i],header=T,sep='\t',stringsAsFactors=F)
	
	# remove duplicate gene entries (i.e. only keep row with highest logFC per gene)
	ipDf <- ipDf[!duplicated(ipDf$gene),]
	
	allDf <- rbind(allDf,data.frame(IP=summaryDf$IP[i],
		ipDf[,c('gene','logFC','FDR','significant')]))
}


# ----------------------------------------------------------------------------------------
# (1) Calculate pairwise similarity metrics for each pair of IPs

# function to calculate Jaccard index of detected proteins or significant interactors
# and overlap enrichment of significant interactors
getJaccardIndex <- function(IP1,IP2) {
	IP1_bait <- subset(summaryDf,IP==IP1)$Bait
	IP2_bait <- subset(summaryDf,IP==IP2)$Bait
	
	# overlap of detected proteins
	IP1_det <- setdiff(unique(subset(allDf,IP==IP1)$gene),IP1_bait)
	IP2_det <- setdiff(unique(subset(allDf,IP==IP2)$gene),IP2_bait)
	overlap_det <- intersect(IP1_det,IP2_det)
	union_det <- union(IP1_det,IP2_det)
	jaccard_det <- length(overlap_det)/length(union_det)

	# overlap of significant interactors
	IP1_ints <- setdiff(unique(subset(allDf,IP==IP1 & significant)$gene),IP1_bait)
	IP2_ints <- setdiff(unique(subset(allDf,IP==IP2 & significant)$gene),IP2_bait)
	overlap_ints <- intersect(IP1_ints,IP2_ints)
	union_ints <- union(IP1_ints,IP2_ints)
	jaccard_ints <- length(overlap_ints)/length(union_ints)

	# overlap of significant interactors (restricted to overlap_det)
	IP1_ints_inter <- intersect(IP1_ints,overlap_det)
	IP2_ints_inter <- intersect(IP2_ints,overlap_det)
	overlap_ints_inter <- intersect(overlap_ints,overlap_det)
	union_ints_inter <- intersect(union_ints,overlap_det)
	jaccard_ints_inter <- length(overlap_ints_inter)/length(union_ints_inter)

	# overlap enrichment (overlap_det as population)
	q <- length(overlap_ints_inter)
	m <- length(IP1_ints_inter)
	n <- length(setdiff(overlap_det,IP1_ints_inter))
	k <- length(IP2_ints_inter)
	p_ints_inter <- phyper(q-1,m,n,k,lower.tail=F)

	# store all stats in data frame
	df <- data.frame(IP1_Det=length(IP1_det),IP2_Det=length(IP2_det),
		Overlap_Det=length(overlap_det),Union_Det=length(union_det),
		Jaccard_Det=jaccard_det,

		IP1_Ints=length(IP1_ints),IP2_Ints=length(IP2_ints),
		Overlap_Ints=length(overlap_ints),Union_Ints=length(union_ints),
		Jaccard_Ints=jaccard_ints,
		
		IP1_Ints_inter=length(IP1_ints_inter),IP2_Ints_inter=length(IP2_ints_inter),
		Overlap_Ints_inter=length(overlap_ints_inter),
		Union_Ints_inter=length(union_ints_inter),
		Jaccard_Ints_inter=jaccard_ints_inter,OverlapP_inter=p_ints_inter,
		Overlap_Ints_List=paste(overlap_ints_inter,collapse=','))
	
	return(df)
}


# function to calculate log2 FC or signed -log10 FDR correlation of detected proteins
getFcFdrCorr <- function(IP1,IP2) {
	IP1_bait <- subset(summaryDf,IP==IP1)$Bait
	IP2_bait <- subset(summaryDf,IP==IP2)$Bait
	IP1_df <- subset(allDf,IP==IP1 & gene != IP1_bait)
	IP2_df <- subset(allDf,IP==IP2 & gene != IP2_bait)
	inter_df <- merge(IP1_df,IP2_df,by='gene')

	fc_corr <- cor(inter_df$logFC.x,inter_df$logFC.y)

	fdr_corr <- cor(sign(inter_df$FDR.x)*-log10(inter_df$FDR.x),
		sign(inter_df$FDR.y)*-log10(inter_df$FDR.y))

	return(data.frame(FC_Corr=fc_corr,FDR_Corr=fdr_corr))
}


# calculate and store similarity metrics for all IP pairs
pairDf <- expand.grid(IP1=summaryDf$IP,IP2=summaryDf$IP)
pairDf <- pairDf %>% rowwise() %>% 
	mutate(getJaccardIndex(IP1,IP2),getFcFdrCorr(IP1,IP2))


# output similarity metrics (excluding symmetrical entries)
# compare similarity metrics of IPs with same vs. different variables
outDf <- NULL
for (i in 1:(length(summaryDf$IP)-1)) {
	for (j in (i+1):length(summaryDf$IP)) {
	
		outDf <- rbind(outDf,subset(pairDf,IP1==summaryDf$IP[i] & IP2==summaryDf$IP[j]))	
	}
}

# output similarilty metrics
write.table(outDf,'../output/CAD_IpPairwiseStats.txt',sep='\t',quote=F,row.names=F)


# ------------------------------------------------------------------------------
# (2) Generate Venn diagrams of interactors for IPs of the same bait

pdf('../output/CAD_SameBaitIps_VennDiagrams.pdf',width=5,height=5)

# iterate through baits with multiple IPs
multiBaits <- names(table(summaryDf$Bait)[table(summaryDf$Bait)>1])
for (bait in multiBaits) {

	baitIps <- subset(summaryDf,Bait==bait)$IP
	df <- subset(allDf,IP %in% baitIps & gene != bait & significant)[,c('IP','gene')]
	df <- aggregate(.~IP,df,c) # get list of interactors in each IP
	df$IP <- sapply(strsplit(df$IP,paste('.',bait,sep='')),'[[',1) # suffix of IP name

	# format interactor lists for input into ggVennDiagram below
	intsList <- vector(mode='list',length=length(df$IP))
	names(intsList) <- df$IP
	for (i in 1:nrow(df)) {
		intsList[[df$IP[i]]] <- unlist(df$gene[i])
	}

	# plot venn diagram
	p <- ggVennDiagram(intsList,label_alpha=0,set_size=3,label_size=3,
		edge_color='black') +
	
		# expand axis to show full IP names
		scale_x_continuous(expand=expansion(mult=0.2)) +
	
		# fill color based on number in each region
		scale_fill_gradient(low='white',high='orange') +
		
		ggtitle(paste(bait,'interactors')) + 
		theme(plot.title = element_text(face='bold',size=10,hjust = 0.5))
		
	print(p)
}

dev.off()
