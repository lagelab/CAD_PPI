##########################################################################################
## Compare expression of interactors, non-interactors, and other proteins
## in whole proteome data from human heart cell types (Doll et al. Nat Commun 2017)
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

# Code tested with R version 4.2.2
library(readxl) # v1.4.1
library(tidyverse) # v1.3.2
library(ggpubr) # v0.5.0

# ----------------------------------------------------------------------------------------
# read in whole proteome data from Supp Data 7 of Doll et al.
# protein expression measured in cardiac fibroblasts (CF), endothelial cells (EC),
# smooth muscle cells (SMC), and adipose fibroblasts (AF)

wpDf <- read_excel('../data/Doll2017_HeartProteome_SuppData7.xlsx',
	sheet='All cell type data')
wpDf <- wpDf[,c('AF','EC','SMC','CF','T: Gene names')]
names(wpDf)[5] <- 'Gene'
wpDf$Gene <- sapply(strsplit(wpDf$Gene,';'),'[[',1)


# ----------------------------------------------------------------------------------------
# read in ints and non-ints in 6 COMBINED PPI networks:
# EC, SMC, Union, EC only, SMC only, Intersect
# extract expression values from proteome data above

netDf <- read.table('../data/CAD_Network-Cell-Bait_Mapping.txt',
	header=T,sep='\t',stringsAsFactors=F)
netDf <- subset(netDf,grepl('COMBINED',Network))
netDf$CellType[netDf$CellType=='EC-SMC'] <- 'Union'
netDf$CellType[netDf$CellType=='EC-SMC-intersect'] <- 'Intersect'
netDf$CellType <- gsub('-',' ',netDf$CellType)

intDf <- read.table('../data/CAD_MasterInteractorTable.txt',
	header=T,sep='\t',stringsAsFactors=F)
intDf <- subset(intDf,intDf$ListName %in% netDf$Network)


expDf <- NULL
# iterate through each network
for (i in 1:nrow(netDf)) {

	# define Int vs. NonInt
	df <- subset(intDf,intDf$ListName==netDf$Network[i])
	df$Type <- ifelse(df$Interactor,'Int','NonInt')
	
	# restrict to gene names found in whole proteome data
	df <- subset(df,GeneName %in% wpDf$Gene)

	# network baits, to be excluded below
	baits <- unlist(strsplit(netDf$Baits[i],','))
	
	for (ct in c('AF','EC','SMC','CF')) {
		# store expression of Int and NonInt
		expVals <- wpDf[match(df$GeneName,wpDf$Gene),ct] # take first match
		expDf <- rbind(expDf,data.frame(Network=netDf$CellType[i],Proteome=ct,
			GeneName=df$GeneName,Type=df$Type,Exp=pull(expVals)))
			
		# store expression of other proteins in proteome data
		otherDf <- subset(wpDf,(!Gene %in% df$GeneName) & (!Gene %in% baits))
		expDf <- rbind(expDf,data.frame(Network=netDf$CellType[i],Proteome=ct,
			GeneName=otherDf$Gene,Type='Other',Exp=pull(otherDf[,ct])))
	}
}


# ----------------------------------------------------------------------------------------
# plot and compare expression of Int, NonInt, and Other proteins

expDf$Network <-factor(expDf$Network,levels=unique(expDf$Network))
expDf$Proteome <-factor(expDf$Proteome,levels=unique(expDf$Proteome))
expDf$Type <- factor(expDf$Type,levels=c('Int','NonInt','Other'))

# set Wilcoxon test parameters
varPairs <- combn(levels(expDf$Type),2)
testList <- NULL
for (col in 1:ncol(varPairs)) { testList[[col]] <- varPairs[,col] }
pCutoffs <- c(0,0.001,0.05,1)
pSymbols <- c('**','*','ns')


pdf('../output/CAD_HeartProteome_ViolinPlots.pdf',height=8,width=6)

ggviolin(expDf,x='Type',y='Exp',color='Type',add='boxplot') +
facet_grid(Network~Proteome) +

# pairwise Wilcoxon tests
stat_compare_means(comparisons=testList,method='wilcox.test',
	symnum.args=list(cutpoints=pCutoffs,symbols=pSymbols),hide.ns=T,size=3,vjust=0.5) +

scale_color_brewer(palette='Dark2') +
xlab('Protein type') + ylab('Relative expression') +
theme_bw() + theme(legend.position='none')

dev.off()
