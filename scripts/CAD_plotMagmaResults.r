##########################################################################################
## Generate plots for CAD PPI MAGMA results
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

# Code tested with R version 4.2.2
library(ggplot2) # v3.4.0
library(tidyr) # v1.3.0

# ----------------------------------------------------------------------------------------
# read in MAGMA results
df <- read.table('../data/CAD_MagmaResults.txt',header=T,stringsAsFactors=F)

# set test plotting order
df$TEST <- factor(df$TEST,levels=c('Global','Conditional'))

# set bait/network plotting order
df$Bait <- sapply(strsplit(df$VARIABLE,'_'),'[',1)
baitOrder <- c('COMBINED',sort(unique(df$Bait[df$Bait!='COMBINED'])))
df$Bait <- factor(df$Bait, levels=rev(baitOrder))

# set cell type plotting order
df$Cell <- sapply(strsplit(df$VARIABLE,'_'),'[',2)
df$Cell[df$Cell=='EC-SMC'] <- 'Union'
df$Cell[df$Cell=='EC-only'] <- 'EC only'
df$Cell[df$Cell=='SMC-only'] <- 'SMC only'
df$Cell[df$Cell=='EC-SMC-intersect'] <- 'Intersect'
cellOrder <- c('EC','Union','SMC','EC only','Intersect','SMC only')
df$Cell <- factor(df$Cell,levels=cellOrder)

# set trait plotting order
df$TRAIT <- gsub('[A-Za-z]+-','',df$TRAIT)
traitOrder <- c('CAD','AA','DA','AS','AIS','LAS','CES','SVS','height')
df$TRAIT <- factor(df$TRAIT, levels=traitOrder)

# set significance marker text for heat map
df$Sig <- ifelse(df$P < 0.05,'*','')
df$Sig[df$P < 0.05/29] <- '**' # adjusting for 29 networks (EC, SMC, Union networks) 


# ----------------------------------------------------------------------------------------
# heat maps without height, EC or SMC global results only
pdf('../output/CAD_MagmaResults_GlobalSingleHeatMaps.pdf',height=3,width=3.5)

# EC
tempDf <- subset(df,TEST=='Global' & TRAIT!='height' & Cell=='EC')

ggplot(tempDf,aes(x=TRAIT,y=Bait,fill=-log10(P))) +
geom_tile() + geom_text(aes(label=Sig)) +
scale_fill_gradient(name=expression(paste(-log[10],"(P-value)",sep="")),
	low='white', high='red',na.value='grey') +
xlab('GWAS phenotype') + ylab('Index protein') + ggtitle('HAEC') +
theme_bw() + 
theme(legend.position="bottom",legend.box="vertical",
	legend.margin=margin(),legend.box.margin=margin(-5,10,-5,-10),
	legend.key.size=unit(0.9,"line"),legend.title=element_text(size=9),
	legend.text=element_text(size=8),axis.text=element_text(size=9))

# SMC
tempDf <- subset(df,TEST=='Global' & TRAIT!='height' & Cell=='SMC')

ggplot(tempDf,aes(x=TRAIT,y=Bait,fill=-log10(P))) +
geom_tile() + geom_text(aes(label=Sig)) +
scale_fill_gradient(name=expression(paste(-log[10],"(P-value)",sep="")),
	low='white', high='red',na.value='grey') +
xlab('GWAS phenotype') + ylab('Index protein') + ggtitle('HCASMC') + 
theme_bw() + 
theme(legend.position="bottom",legend.box="vertical",
	legend.margin=margin(),legend.box.margin=margin(-5,10,-5,-10),
	legend.key.size=unit(0.9,"line"),legend.title=element_text(size=9),
	legend.text=element_text(size=8),axis.text=element_text(size=9))

dev.off()


# ----------------------------------------------------------------------------------------
# fill in NA for empty cells (to be plotted as grey in heat map)
df <- complete(df,TEST,Cell,TRAIT,Bait) 

# heat maps with cell type as facets
pdf('../output/CAD_MagmaResults_HeatMapFacet.pdf',height=5,width=6.5)

# global
ggplot(subset(df,TEST=='Global'),aes(x=TRAIT,y=Bait,fill=-log10(P))) +
facet_wrap(~Cell) + geom_tile() + geom_text(aes(label=Sig)) +
scale_fill_gradient(name=expression(paste(-log[10],"(P-value)",sep="")),
	low='white', high='red',na.value='grey') +
xlab('GWAS phenotype') + ylab('Index protein') + ggtitle('Global') +
theme_bw() + 
theme(legend.position="bottom",legend.box="vertical",
	legend.margin=margin(),legend.box.margin=margin(-5,10,-5,-10),
	legend.key.size=unit(0.9,"line"),legend.title=element_text(size=9),
	legend.text=element_text(size=8),
	axis.text.x=element_text(size=9,angle=45,hjust=0.9),axis.text.y=element_text(size=9))

# conditional
ggplot(subset(df,TEST=='Conditional'),aes(x=TRAIT,y=Bait,fill=-log10(P))) +
facet_wrap(~Cell) + geom_tile() + geom_text(aes(label=Sig)) +
scale_fill_gradient(name=expression(paste(-log[10],"(P-value)",sep="")),
	low='white', high='red',na.value='grey') +
xlab('GWAS phenotype') + ylab('Index protein') + ggtitle('Conditional') +
theme_bw() + 
theme(legend.position="bottom",legend.box="vertical",
	legend.margin=margin(),legend.box.margin=margin(-5,10,-5,-10),
	legend.key.size=unit(0.9,"line"),legend.title=element_text(size=9),
	legend.text=element_text(size=8),
	axis.text.x=element_text(size=9,angle=45,hjust=0.9),axis.text.y=element_text(size=9))

dev.off()


# ----------------------------------------------------------------------------------------
# heat maps without height, with cell type as facets
pdf('../output/CAD_MagmaResults_HeatMapFacet_NoHeight.pdf',height=5,width=6)

# global
ggplot(subset(df,TEST=='Global' & TRAIT!='height'),aes(x=TRAIT,y=Bait,fill=-log10(P))) +
facet_wrap(~Cell) + geom_tile() + geom_text(aes(label=Sig)) +
scale_fill_gradient(name=expression(paste(-log[10],"(P-value)",sep="")),
	low='white', high='red',na.value='grey') +
xlab('GWAS phenotype') + ylab('Index protein') + ggtitle('Global') +
theme_bw() + 
theme(legend.position="bottom",legend.box="vertical",
	legend.margin=margin(),legend.box.margin=margin(-5,10,-5,-10),
	legend.key.size=unit(0.9,"line"),legend.title=element_text(size=9),
	legend.text=element_text(size=8),
	axis.text.x=element_text(size=9,angle=45,hjust=0.9),axis.text.y=element_text(size=9))

# conditional
ggplot(subset(df,TEST=='Conditional' & TRAIT!='height'),
	aes(x=TRAIT,y=Bait,fill=-log10(P))) +
facet_wrap(~Cell) + geom_tile() + geom_text(aes(label=Sig)) +
scale_fill_gradient(name=expression(paste(-log[10],"(P-value)",sep="")),
	low='white', high='red',na.value='grey') +
xlab('GWAS phenotype') + ylab('Index protein') + ggtitle('Conditional') +
theme_bw() + 
theme(legend.position="bottom",legend.box="vertical",
	legend.margin=margin(),legend.box.margin=margin(-5,10,-5,-10),
	legend.key.size=unit(0.9,"line"),legend.title=element_text(size=9),
	legend.text=element_text(size=8),
	axis.text.x=element_text(size=9,angle=45,hjust=0.9),axis.text.y=element_text(size=9))

dev.off()
