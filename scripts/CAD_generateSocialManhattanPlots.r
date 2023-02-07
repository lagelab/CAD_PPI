##########################################################################################
## Generate CAD PPI social Manhattan plots
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

# Code tested with R version 4.2.2
library(readxl) # 1.4.1
library(ggplot2) # 3.4.0
library(dplyr) # 1.1.0
library(ggrepel) # 0.9.3
library(ggnewscale) # 0.4.8

# ----------------------------------------------------------------------------------------
# read in genes to be plotted in social Manhattan plots
ecDf <- read_excel("../data/CAD_SocialManhattanData.xlsx",sheet='EC')
smcDf <- read_excel("../data/CAD_SocialManhattanData.xlsx",sheet='SMC')

# merge EC and SMC interactions, add CellType column
jointDf <- data.frame(intersect(ecDf,smcDf),CellType='EC-SMC')
jointDf <- rbind(jointDf,data.frame(setdiff(ecDf,smcDf),CellType='EC'))
jointDf <- rbind(jointDf,data.frame(setdiff(smcDf,ecDf),CellType='SMC'))
dim(jointDf)

# list of the 3 data frams
dfList <- list('EC'=ecDf,'SMC'=smcDf,'Joint'=jointDf)

for (d in names(dfList)) {

	df <- as.data.frame(dfList[[d]])

	# set chr 1-22 coordinates on x-axis
	df$Coord <- NA
	s <- 0
	
	# GRCh37 chr 1-22 lengths
	nbp <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,
		146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,
		90354753,81195210,78077248,59128983,63025520,48129895,51304566)
		
	chrTicks <- NULL
	for (i in 1:22) {
		df[df$Chr==i,"Coord"] <- df[df$Chr==i,"Pos"] + s
		chrTicks <- c(chrTicks,floor(nbp[i]/2) + s)
		s <- s + nbp[i]
	}

	# add bait coordinates (for plotting social links between genes and baits)
	df$Bait_P <- NA
	df$Bait_Coord <- NA
	for (bait in unique(df$Bait)) {

		df[df$Bait==bait,"Bait_P"] <- unique(df$P[df$Gene==bait])
		df[df$Bait==bait,"Bait_Coord"] <- unique(df$Coord[df$Gene==bait])
		# unique() to get rid of duplicate values (if bait shows up as interactor)

	}
	# don't plot self link to bait itself
	df[df$Gene==df$Bait,c('Bait_P','Bait_Coord')]<- NA
	
	dfList[[d]] <- df
}


pdf("../output/CAD_SocialManhattanPlots.pdf",width=6.5,height=3.5)

# ----------------------------------------------------------------------------------------
# JOINT (EC+SMC) PLOT

df <- dfList[['Joint']]

# node color
df$Group <- ifelse(df$Gene %in% unique(df$Bait),'Index','Locus')
# WB = interactors replicated by western blotting
df$Group[df$Gene %in% c('IGF2BP1','MAP4','ATXN2','TNS1','FNDC3B')] <- 'WB'
df$Group <- factor(df$Group,levels=c('Locus','WB','Index'))
df <- df[order(df$Group),] # to plot index gene data points last

# edge color
df$CellType <- factor(df$CellType,levels=c('EC','SMC','EC-SMC'))

# unique genes that show up in df (which contains unique PAIRS of genes)
uniqDf <- unique(df[,c("Gene","Coord","P","Group")])
table(uniqDf$Group)
# 11 index genes + 61 locus genes

set.seed(123) # fix ggrepel label coordinates

ggplot(uniqDf) + 

# social links
geom_segment(df,mapping=aes(x=Coord,y=-log10(P),xend=Bait_Coord,yend=-log10(Bait_P),
	color=CellType),size=0.3,alpha=0.6,show.legend=F) +

scale_color_manual(labels=c('EC only','SMC only','Both'),
	values=c('dodgerblue','coral','black')) +

# genes
ggnewscale::new_scale_color() +
geom_point(uniqDf,mapping=aes(x=Coord,y=-log10(P),color=Group),size=1,show.legend=F) +

geom_text_repel(uniqDf,mapping=aes(x=Coord,y=-log10(P),label=Gene,color=Group),
	size=ifelse(uniqDf$Group=='Index',2,1.6),
	fontface='bold',segment.size=0.1,segment.alpha=0.5,
	max.overlaps=100,show.legend=FALSE) +

scale_color_manual(labels=c('Locus','WB','Index'),values=c('black','magenta','red')) +

scale_x_continuous(limits=c(0,max(df$Coord)),label=as.character(1:22),breaks=chrTicks) +
scale_y_continuous(expand=c(0,0), limits=c(0,31)) +
xlab("Chromosomal position") + ylab(expression(paste(-log[10],"(GWAS P-value)",sep=""))) +
ggtitle('HAEC + HCASMC') + theme_classic() +
theme(axis.title=element_text(size=8),axis.text=element_text(size=7))


# ----------------------------------------------------------------------------------------
# JOINT (EC+SMC) PLOT showing only EDN1 or PHACTR1 interactions

df <- subset(dfList[['Joint']],Bait=='EDN1' | Bait=='PHACTR1')

# node color
df$Group <- 'Locus'
df$Group[df$Gene=='EDN1'] <- 'EDN1'
df$Group[df$Gene=='PHACTR1'] <- 'PHACTR1'
df$Group <- factor(df$Group,levels=c('Locus','PHACTR1','EDN1'))
df <- df[order(df$Group),] # to plot index gene data points last

# unique genes that show up in df (which contains unique PAIRS of genes)
uniqDf <- unique(df[,c("Gene","Coord","P","Group")])

set.seed(123) # fix ggrepel label coordinates

ggplot(uniqDf) + 

# social links
geom_segment(df,mapping=aes(x=Coord,y=-log10(P),xend=Bait_Coord,yend=-log10(Bait_P)),
	color=ifelse(df$Bait=='EDN1','red','blue'),size=0.5,alpha=0.4,show.legend=F) +

# genes
geom_point(uniqDf,mapping=aes(x=Coord,y=-log10(P),color=Group),size=1,show.legend=F) +

geom_text_repel(uniqDf,mapping=aes(x=Coord,y=-log10(P),label=Gene,color=Group),
	size=2,fontface='bold',segment.size=0.1,segment.alpha=0.5,
	max.overlaps=100,show.legend=FALSE) +

scale_color_manual(labels=c('Locus','PHACTR1','EDN1'),values=c('black','blue','red')) +

scale_x_continuous(limits=c(0,max(df$Coord)),label=as.character(1:22),breaks=chrTicks) +
scale_y_continuous(expand=c(0,0), limits=c(0,31)) +
xlab("Chromosomal position") + ylab(expression(paste(-log[10],"(GWAS P-value)",sep=""))) +
ggtitle('HAEC + HCASMC, EDN1 vs. PHACTR1 interactions') + theme_classic() +
theme(axis.title=element_text(size=8),axis.text=element_text(size=7))


# ----------------------------------------------------------------------------------------
# EC PLOT

df <- dfList[['EC']]

# node color
df$Group <- ifelse(df$Gene %in% unique(df$Bait),'Index','Locus')
df$Group <- factor(df$Group,levels=c('Locus','Index'))
df <- df[order(df$Group),] # to plot index gene data points last

# edge color (interactions replicated by western blotting)
df$Validated <- FALSE
df[rownames(subset(df,(Gene=='IGF2BP1' & Bait=='BCAS3')|(Gene=='MAP4' & Bait=='EDN1')|
	(Gene=='MAP4' & Bait=='FN1')|(Gene=='IGF2BP1' & Bait=='PLPP3')|
	(Gene=='MAP4' & Bait=='PLPP3'))),'Validated'] <- TRUE
df <- df[order(df$Validated),]

# unique genes that show up in df (which contains unique PAIRS of genes)
uniqDf <- unique(df[,c("Gene","Coord","P","Group")])
table(uniqDf$Group)

set.seed(123) # fix ggrepel label coordinates

ggplot(uniqDf) + 

# social links
geom_segment(df,mapping=aes(x=Coord,y=-log10(P),xend=Bait_Coord,yend=-log10(Bait_P),
	color=Validated),size=0.3,alpha=0.6,show.legend=F) +

scale_color_manual(labels=c('FALSE','TRUE'),values=c('grey','blue')) +

# genes
ggnewscale::new_scale_color() +
geom_point(uniqDf,mapping=aes(x=Coord,y=-log10(P),color=Group),size=1,show.legend=F) +

geom_text_repel(uniqDf,mapping=aes(x=Coord,y=-log10(P),label=Gene,color=Group),
	size=ifelse(uniqDf$Group=='Index',2,1.6),
	fontface='bold',segment.size=0.1,segment.alpha=0.5,
	max.overlaps=100,show.legend=FALSE) +

scale_color_manual(labels=c('Locus','Index'),values=c('black','red')) +

scale_x_continuous(limits=c(0,max(df$Coord)),label=as.character(1:22),breaks=chrTicks) +
scale_y_continuous(expand=c(0,0), limits=c(0,31)) +
xlab("Chromosomal position") + ylab(expression(paste(-log[10],"(GWAS P-value)",sep=""))) +
ggtitle('HAEC') + theme_classic() +
theme(axis.title=element_text(size=8),axis.text=element_text(size=7))


# ----------------------------------------------------------------------------------------
# SMC PLOT

df <- dfList[['SMC']]

# node color
df$Group <- ifelse(df$Gene %in% unique(df$Bait),'Index','Locus')
df$Group <- factor(df$Group,levels=c('Locus','Index'))
df <- df[order(df$Group),] # to plot index gene data points last

# edge color (interactions replicated by western blotting)
df$Validated <- FALSE
df[rownames(subset(df,(Gene=='ATXN2' & Bait=='ADAMTS7') |
	(Gene=='IGF2BP1' & Bait=='ADAMTS7') | (Gene=='MAP4' & Bait=='EDN1') |
	(Gene=='FNDC3B' & Bait=='EDNRA') | (Gene=='TNS1' & Bait=='EDNRA') |
	(Gene=='IGF2BP1' & Bait=='FN1') | (Gene=='MAP4' & Bait=='FN1') |
	(Gene=='FNDC3B' & Bait=='JCAD') | (Gene=='IGF2BP1' & Bait=='JCAD') |
	(Gene=='MAP4' & Bait=='JCAD') | (Gene=='TNS1' & Bait=='JCAD') |
	(Gene=='MAP4' & Bait=='PLPP3') | (Gene=='TNS1' & Bait=='PLPP3'))),'Validated'] <- TRUE
df <- df[order(df$Validated),]

# unique genes that show up in df (which contains unique PAIRS of genes)
uniqDf <- unique(df[,c("Gene","Coord","P","Group")])
table(uniqDf$Group)

set.seed(42) # fix ggrepel label coordinates

ggplot(uniqDf) + 

# social links
geom_segment(df,mapping=aes(x=Coord,y=-log10(P),xend=Bait_Coord,yend=-log10(Bait_P),
	color=Validated),size=0.3,alpha=0.6,show.legend=F) +

scale_color_manual(labels=c('FALSE','TRUE'),values=c('grey','blue')) +

# genes
ggnewscale::new_scale_color() +
geom_point(uniqDf,mapping=aes(x=Coord,y=-log10(P),color=Group),size=1,show.legend=F) +

geom_text_repel(uniqDf,mapping=aes(x=Coord,y=-log10(P),label=Gene,color=Group),
	size=ifelse(uniqDf$Group=='Index',2,1.6),
	fontface='bold',segment.size=0.1,segment.alpha=0.5,
	max.overlaps=100,show.legend=FALSE) +

scale_color_manual(labels=c('Locus','Index'),values=c('black','red')) +

scale_x_continuous(limits=c(0,max(df$Coord)),label=as.character(1:22),breaks=chrTicks) +
scale_y_continuous(expand=c(0,0), limits=c(0,31)) +
xlab("Chromosomal position") + ylab(expression(paste(-log[10],"(GWAS P-value)",sep=""))) +
ggtitle('HCASMC') + theme_classic() +
theme(axis.title=element_text(size=8),axis.text=element_text(size=7))
	
dev.off()
