##########################################################################################
## UGER job submission script
## run global MAGMA analysis for each CAD PPI network (ints vs. genomic background)
## Step 1 of 2: run gene annotation step for each network (same across all traits)
##
## Author: Yu-Han Hsu
##########################################################################################

#$ -cwd
#$ -N uger.magma.global.annot
#$ -l h_vmem=4g
#$ -l h_rt=04:00:00
#$ -t 1-57
#$ -tc 57


listName=$(awk "NR==$SGE_TASK_ID+1 {print \$1}" CAD_Network-Cell-Bait_Mapping.txt)
echo $listName

# create gene.loc file for each network
mkdir -p temp

awk "NR==$SGE_TASK_ID+1 {print \$3}" CAD_Network-Cell-Bait_Mapping.txt | tr ',' '\n' > \
temp/${listName}.baits.txt

grep -vwf temp/${listName}.baits.txt Ensembl_BioMart_GRCh37.gene.loc > \
temp/${listName}.GLOBAL.gene.loc


# directory to store MAGMA output files
mkdir -p magma_output

# MAGMA directory
magmaDir=Software/MAGMA

# Annotation Step (SNP to gene mapping)
${magmaDir}/magma_v1.09_static/magma \
--annotate window=50 \
filter=${magmaDir}/ReferenceData/g1000_eur/g1000_eur.noMHC.bim \
--snp-loc ${magmaDir}/ReferenceData/g1000_eur/g1000_eur.bim \
--gene-loc temp/${listName}.GLOBAL.gene.loc \
--out magma_output/${listName}.GLOBAL


# delete temp files
rm temp/${listName}.baits.txt temp/${listName}.GLOBAL.gene.loc
