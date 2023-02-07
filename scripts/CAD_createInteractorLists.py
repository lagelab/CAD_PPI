##########################################################################################
## Define lists of interactors and non-interactors for EC, SMC, and Union PPI networks
##
## Author: Yu-Han Hsu
##########################################################################################

# Code tested using Python 3.8.9
import pandas as pd # v1.4.2 (check using: pd. __version__)


# output file path and header
outFile = open('../output/CAD_EC_SMC_Union_InteractorTable.txt','w')
outFile.write('ListName\tBait\tCellType\tGeneName\tInteractor\tChr\tStart\tEnd\n')

# gene chrmosomal position dict
posDict = {} # gene name as key, [chr, start, end] (hg19) as value
with open('../data/CAD_gene_hg19_positions.tsv','r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.strip().split('\t')
		posDict[fields[0]] = fields[1:4]

# list of 20 IPs that passed QC
ipFiles = []
baits = []
cells = [] # EC or SMC
with open('../data/CAD_IpSummary.txt','r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.strip().split('\t')
		ipFiles.append(fields[8])
		baits.append(fields[1])
		cells.append(fields[2])

listNames = [x+'_'+y for x,y in zip(baits, cells)]

# ----------------------------------------------------------------------------------------
# create interactor lists for individual BAIT_CELLTYPE networks
for listName in list(set(listNames)):

	print('##########')
	print(listName)
	
	bait = listName.split('_')[0]
	cell = listName.split('_')[1]
    
    # get all relevant IP files
	baitFiles = [f for f in ipFiles if (bait in f and cell in f)]
	intSet = set()
	nonIntSet = set()

	# iterate through each IP file
	for f in baitFiles:
	
		print('input: ' + f)
		ipTable = pd.read_csv(f,sep="\t")

		# confirm bait name is correct/found in data file
		if (ipTable['gene']==bait).sum() == 0: print('CHECK BAIT NAME')

		# add interactors (logFC > 0 & FDR ≤ 0.1)
		intInds = ((ipTable['gene']!=bait) & 
			(ipTable['logFC'] >= 0) & (ipTable['FDR'] <= 0.1))
		intSet = intSet | set(ipTable[intInds]['gene'])

		# add non-interactors
		nonIntInds = (ipTable['gene']!=bait) & ~intInds
		nonIntSet = nonIntSet | set(ipTable[nonIntInds]['gene'])
	
	# for gene names in both sets, keep in intSet only
	nonIntSet = nonIntSet - intSet

	print('# interactors: %d' % len(intSet))
	print('# non-interactors: %d' % len(nonIntSet))

	# output interactor table with chromosomal position
	for x in intSet:
		if x in posDict:
			outInfo = [listName,bait,cell,x,'TRUE'] + posDict[x]
			outFile.write('\t'.join(outInfo)+'\n')
		else:
			print('Interactor not in posDict: ' + x)

	for x in nonIntSet:
		if x in posDict:
			outInfo = [listName,bait,cell,x,'FALSE'] + posDict[x]
			outFile.write('\t'.join(outInfo)+'\n')
		else:
			print('Non-interactor not in posDict: ' + x)


# ----------------------------------------------------------------------------------------
# create interactor lists for cell-type-combined networks of individual baits
cell = 'EC-SMC' # i.e. "Union"
for bait in list(set(baits)):

	# only do this for baits with data in both cell types
	if (bait + '_EC' in listNames) and (bait + '_SMC' in listNames): 
		listName = bait + '_' + cell
		
		print('##########')
		print(listName)

		# get all relevant IP files
		baitFiles = [f for f in ipFiles if bait in f]
		intSet = set()
		nonIntSet = set()
		
		# iterate through each IP file
		for f in baitFiles:	

			print('input: ' + f)
			ipTable = pd.read_csv(f,sep="\t")

			# add interactors (logFC > 0 & FDR ≤ 0.1)
			intInds = ((ipTable['gene']!=bait) & 
				(ipTable['logFC'] >= 0) & (ipTable['FDR'] <= 0.1))
			intSet = intSet | set(ipTable[intInds]['gene'])

			# add non-interactors
			nonIntInds = (ipTable['gene']!=bait) & ~intInds
			nonIntSet = nonIntSet | set(ipTable[nonIntInds]['gene'])
		
		# for gene names in both sets, keep in intSet only	
		nonIntSet = nonIntSet - intSet
		
		print('# interactors: %d' % len(intSet))
		print('# non-interactors: %d' % len(nonIntSet))

		# output interactor table with chromosomal position
		for x in intSet:
			if x in posDict:
				outInfo = [listName,bait,cell,x,'TRUE'] + posDict[x]
				outFile.write('\t'.join(outInfo)+'\n')

		for x in nonIntSet:
			if x in posDict:
				outInfo = [listName,bait,cell,x,'FALSE'] + posDict[x]
				outFile.write('\t'.join(outInfo)+'\n')


# ----------------------------------------------------------------------------------------
# create interactor lists for bait-combined network in EC or SMC
bait = 'COMBINED'
for cell in list(set(cells)):
	listName = bait + '_' + cell

	print('##########')
	print(listName)

	# get all relevant IP files
	cellFiles = [f for f in ipFiles if cell in f]
	intSet = set()
	nonIntSet = set()
	
	# iterate through each IP file
	for f in cellFiles:
		
		print('input: ' + f)
		ipTable = pd.read_csv(f,sep="\t")

		# add interactors (logFC > 0 & FDR ≤ 0.1)
		intInds = ((~ipTable['gene'].isin(baits)) & 
			(ipTable['logFC'] >= 0) & (ipTable['FDR'] <= 0.1))
		intSet = intSet | set(ipTable[intInds]['gene'])

		# add non-interactors
		nonIntInds = (~ipTable['gene'].isin(baits)) & ~intInds
		nonIntSet = nonIntSet | set(ipTable[nonIntInds]['gene'])
	
	# for gene names in both sets, keep in intSet only	
	nonIntSet = nonIntSet - intSet
	
	print('# interactors: %d' % len(intSet))
	print('# non-interactors: %d' % len(nonIntSet))

	# output interactor table with chromosomal position
	for x in intSet:
		if x in posDict:
			outInfo = [listName,bait,cell,x,'TRUE'] + posDict[x]
			outFile.write('\t'.join(outInfo)+'\n')

	for x in nonIntSet:
		if x in posDict:
			outInfo = [listName,bait,cell,x,'FALSE'] + posDict[x]
			outFile.write('\t'.join(outInfo)+'\n')


# ----------------------------------------------------------------------------------------
# create interactor lists for bait and cell-type-combined network
bait = 'COMBINED'
cell = 'EC-SMC' # i.e. "Union"
listName = bait + '_' + cell

print('##########')
print(listName)

intSet = set()
nonIntSet = set()
	
for f in ipFiles:

	print('input: ' + f)
	ipTable = pd.read_csv(f,sep="\t")

	# add interactors (logFC > 0 & FDR ≤ 0.1)
	intInds = ((~ipTable['gene'].isin(baits)) & 
		(ipTable['logFC'] >= 0) & (ipTable['FDR'] <= 0.1))
	intSet = intSet | set(ipTable[intInds]['gene'])

	# add non-interactors
	nonIntInds = (~ipTable['gene'].isin(baits)) & ~intInds
	nonIntSet = nonIntSet | set(ipTable[nonIntInds]['gene'])		

# for gene names in both sets, keep in intSet only	
nonIntSet = nonIntSet - intSet

print('# interactors: %d' % len(intSet))
print('# non-interactors: %d' % len(nonIntSet))


# output interactor table with chromosomal position
for x in intSet:
	if x in posDict:
		outInfo = [listName,bait,cell,x,'TRUE'] + posDict[x]
		outFile.write('\t'.join(outInfo)+'\n')

for x in nonIntSet:
	if x in posDict:
		outInfo = [listName,bait,cell,x,'FALSE'] + posDict[x]
		outFile.write('\t'.join(outInfo)+'\n')

outFile.close()

print('SCRIPT COMPLETED')
