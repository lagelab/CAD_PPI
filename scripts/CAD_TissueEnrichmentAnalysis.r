##########################################################################################
## Tissue enrichment analysis of CAD PPI networks
## using tissue-specific genes from GTEx (based on RNA or protein expression)
## comparisons: 
## 	(1) interactors vs. all GTEx genes (global)
##	(2) interactors vs. non-interactors (conditional)
## 	(3) non-interactors vs. all GTEx genes (as additional control analysis)
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

# Code tested with R version 4.2.2
library(genoppi) # v1.0.13 (for loading GTEx tissue specificity data)
library(ggplot2) # v3.4.0
library(RColorBrewer) # v1.1.3

# ----------------------------------------------------------------------------------------
# read in ints and non-ints in 6 COMBINED PPI networks:
# EC, SMC, Union, EC only, SMC only, Intersect

netDf <- read.table('../data/CAD_Network-Cell-Bait_Mapping.txt',
	header=T,sep='\t',stringsAsFactors=F)
netDf <- subset(netDf,grepl('COMBINED',Network))
netDf$CellType[netDf$CellType=='EC-SMC'] <- 'Union'
netDf$CellType[netDf$CellType=='EC-SMC-intersect'] <- 'Intersect'
netDf$CellType <- gsub('-',' ',netDf$CellType)

intDf <- read.table('../data/CAD_MasterInteractorTable.txt',
	header=T,sep='\t',stringsAsFactors=F)
intDf <- subset(intDf,intDf$ListName %in% netDf$Network)


# ----------------------------------------------------------------------------------------
# tissue enrichment analysis using tissue specificity data in genoppi:
#	gtex_rna (Finucane et al. Nat Genet 2018)
# 	gtex_protein (Jiang et al. Cell 2020)

otherDf <- rbind(data.frame(data_source='RNA',gtex_rna),
	data.frame(data_source='protein',gtex_protein))

# tissue category mapping info (from Finucane et al.)
catDf <- read.table('../data/GTEx_TissueCategoryMapping.txt',
	header=T,sep='\t',stringsAsFactors=F)
catDf$DataSource <- sapply(strsplit(catDf$DataSource,'_'),'[[',2)

# run 3 comparisons:
# (1) ints vs. all, (2) ints vs. non-ints, (3) non-ints vs. all
params <- data.frame(NetType=c('Interactors','Interactors','Non-interactors'),
	Background=c('Global','Conditional','Global'))

outDf <- NULL
# iterate through networks
for (net in netDf$Network) {
	print(net)

	# interactors and non-interactors in network
	ints <- subset(intDf,ListName==net & Interactor)$GeneName
	nonInts <- subset(intDf,ListName==net & !Interactor)$GeneName

	# network baits, to be excluded from population below
	baits <- unlist(strsplit(subset(netDf,Network==net)$Baits,','))
			
			
	# iterate through data sources
	for (ds in unique(otherDf$data_source)) {
		
		# all genes included in the data source
		allGenes <- unique(subset(otherDf,data_source==ds)$gene)
		
		
		# iterate through tissues
		for (t in unique(subset(otherDf,data_source==ds)$tissue)) {

			# tissue-specific genes
			tissueGenes <- unique(subset(otherDf,
				data_source==ds & tissue==t & significant)$gene)
			
			# iterate through the 3 comparisons
			for (i in 1:nrow(params)) {
				netType <- params$NetType[i]
				background <- params$Background[i]
				
				# define global or conditional background
				if (background=='Global') { # Global: all genes in data source
					popGenes <- setdiff(allGenes,baits)
				} else { # Conditional: ints + non-ints in data source
					popGenes <- intersect(union(ints,nonInts),allGenes)
				}
				
				# test for enrichment of ints or non-ints
				if (netType=='Interactors') {
					xGenes <- intersect(ints,popGenes)
				} else {
					xGenes <- intersect(nonInts,popGenes)
				}
				
				yGenes <- intersect(tissueGenes,popGenes)
				overlapGenes <- intersect(xGenes,yGenes)
			
				# one-tailed hypergeometric test
				m <- length(xGenes)
				n <- length(setdiff(popGenes,xGenes))
				k <- length(yGenes)
				q <- length(overlapGenes)
				p <- phyper(q-1,m,n,k,lower.tail=F)
			
				# store enrichment results
				outDf <- rbind(outDf,data.frame(
					Network=subset(netDf,Network==net)$CellType,
					DataSource=ds,
					Tissue=subset(catDf,DataSource==ds & Tissue==t)$TissueName,
					TissueCategory=subset(catDf,
						DataSource==ds & Tissue==t)$TissueCategory,
					NetType=netType,Background=background,
					PopCount=length(popGenes),NetCount=m,TissueCount=k,OverlapCount=q,
					OverlapP=p,OverlapGenes=paste(overlapGenes,collapse=',')))
			}
		}
	}
}

# output results table
write.table(outDf,'../output/CAD_TissueEnrichmentResults.txt',
	sep='\t',row.names=F,quote=F)


# ----------------------------------------------------------------------------------------
# enrichment bar plots, colored/ordered by tissue category

# network plotting order
outDf$Network <- factor(outDf$Network,
	levels=c('EC','Union','SMC','EC only','Intersect','SMC only'))

# tissue plotting order: by category, then by name
outDf <- outDf[order(outDf$TissueCategory,outDf$Tissue),]
outDf$Tissue <- factor(outDf$Tissue,levels=rev(unique(outDf$Tissue)))
		
# assign colors to each tissue category
cats <- unique(outDf$TissueCategory)
colors <- brewer.pal(length(cats),'Set1')
names(colors) <- cats


pdf('../output/CAD_TissueEnrichment_Category_BarPlots.pdf',height=4,width=7)
	
for (i in 1:nrow(params)) {

	for (ds in c('RNA','protein')) {
		df <- subset(outDf,NetType==params$NetType[i] & 
			Background==params$Background[i] & DataSource==ds)

		nSets <- length(unique(df$Tissue)) # number of tissues
	
		p <- ggplot(df,aes(x=-log10(OverlapP),y=Tissue,fill=TissueCategory)) + 
			geom_col() + facet_wrap(~Network,ncol=3) +
		
			# line indicating p < 0.05/# of tissues
			geom_vline(xintercept=-log10(0.05/nSets),linetype='dashed',linewidth=0.25) +
		
			scale_fill_manual(values=colors) +
		
			xlab(bquote(-log[10]*"("*italic(.("P"))*"-value)")) +
			ylab('Tissue') +
			ggtitle(paste(ds,': ',params$NetType[i],', ',params$Background[i],sep='')) +
			
			theme_bw() +
			theme(plot.title.position='plot',plot.title=element_text(size=8),
            	legend.title=element_text(size=8),legend.text=element_text(size=7),
				axis.title=element_text(size=8),
				axis.text.y=element_blank(),axis.text.x=element_text(size=7))

		print(p)	
	}
}

dev.off()


# ----------------------------------------------------------------------------------------
# cardiovascular tissue enrichment bar plots, colored/ordered by p-values

pdf('../output/CAD_TissueEnrichment_Cadiovascular_BarPlots.pdf',height=1.2,width=2)

for (i in 1:nrow(params)) {

	for (ds in c('RNA','protein')) {
	
		for (net in unique(outDf$Network)) {
			
			df <- subset(outDf,Network==net & NetType==params$NetType[i] & 
				Background==params$Background[i] & DataSource==ds)
			
			nSets <- length(unique(df$Tissue)) # number of tissues
			df$Sig <- ifelse(df$OverlapP < 0.05,'Nominal','None')
			df$Sig[df$OverlapP < 0.05/nSets] <- 'Bonferroni'

			df <- df[order(df$OverlapP),]
			df$Tissue <- factor(df$Tissue,levels=rev(df$Tissue))

			p <- ggplot(subset(df,TissueCategory=='Cardiovascular'),
				aes(x=-log10(OverlapP),y=Tissue,fill=Sig)) + geom_col() +
		
				# line indicating p < 0.05 and p < 0.05/# of tissues
				geom_vline(xintercept=-log10(0.05),linetype='dashed',linewidth=0.25) +
				geom_vline(xintercept=-log10(0.05/nSets),
					linetype='dashed',linewidth=0.25) +
		
				scale_fill_manual(values=c('None'='black',
					'Nominal'='sienna1','Bonferroni'='red')) +
		
				xlab(bquote(-log[10]*"("*italic(.("P"))*")")) +
				ylab('Tissue') +
				ggtitle(paste(ds,': ',net,' ',params$NetType[i],
					', ',params$Background[i],sep='')) +
			
				theme_bw() +
				theme(plot.title.position='plot',plot.title=element_text(size=6),
					legend.position='none',
					axis.title=element_text(size=6),axis.text=element_text(size=6))

			print(p)
		
		}
	}
}

dev.off()
