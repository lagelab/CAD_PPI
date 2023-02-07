##########################################################################################
## For baits with data in both cell types:
## (1) Define lists of interactors and non-interactors for EC only, SMC only,
## and Intersect PPI networks
## (2) Print int and non-int counts in all networks (including EC, SMC, Union networks)
##
##
## Author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

# Code tested with R version 4.2.2

# read in data for EC, SMC, and Union networks
allDf <- read.table('../data/CAD_EC_SMC_Union_InteractorTable.txt',
	header=T,sep='\t',stringsAsFactors=F)

# only do this for baits with data in both cell types
baits <- c('COMBINED','ARHGEF26','BCAS3','EDN1','FN1','HDAC9','JCAD','PHACTR1','PLPP3')

outDf <- NULL
for (b in baits) {
	print(b)

	ecDf <- subset(allDf,Bait==b & CellType=='EC')
	smcDf <- subset(allDf,Bait==b & CellType=='SMC')
	bothDf <- subset(allDf,Bait==b & CellType=='EC-SMC')
	
	print('EC')
	print(table(ecDf$Interactor))

	print('EC-only')
	eoDf <- subset(ecDf, (Interactor & !(GeneName %in% smcDf$GeneName[smcDf$Interactor]))
		| !Interactor)
	eoDf$ListName <- paste(b,'_EC-only',sep='')
	print(table(eoDf$Interactor))

	print('SMC')
	print(table(smcDf$Interactor))
	
	print('SMC-only')
	soDf <- subset(smcDf, (Interactor & !(GeneName %in% ecDf$GeneName[ecDf$Interactor])) 
		| !Interactor)
	soDf$ListName <- paste(b,'_SMC-only',sep='')
	print(table(soDf$Interactor))

	print('EC-SMC (Union)')
	print(table(bothDf$Interactor))
	
	print('EC-SMC-intersect (Intersect)')
	inDf <- subset(bothDf, (Interactor & (GeneName %in% ecDf$GeneName[ecDf$Interactor]) & 
		(GeneName %in% smcDf$GeneName[smcDf$Interactor])) | !Interactor)
	inDf$ListName <- paste(b,'_EC-SMC-intersect',sep='')
	print(table(inDf$Interactor))

	# store EC only, SMC only, Intersect networks
	outDf <- rbind(outDf,eoDf,soDf,inDf)
}

# output lists of ints and non-ints in EC only, SMC only, Intersect networks
write.table(outDf,'../output/CAD_EConly_SMConly_Intersect_InteractorTable.txt',
	sep='\t',quote=F,row.names=F)
