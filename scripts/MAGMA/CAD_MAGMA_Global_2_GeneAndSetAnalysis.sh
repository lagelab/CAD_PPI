##########################################################################################
## UGER job submission script
## run global MAGMA analysis for each CAD PPI network (ints vs. genomic background)
## Step 2 of 2: run gene and gene set annalysis steps for each network and each trait
## for GWAS traits: CAD, Aorta-AA, Aorta-DA, Stroke-AS, Stroke-AIS, Stroke-LAS,
##                  Stroke-CES, Stroke-SVS, height
##
## Author: Yu-Han Hsu
##########################################################################################

#$ -cwd
#$ -N uger.magma.global.gene
#$ -l h_vmem=4g
#$ -l h_rt=48:00:00
#$ -t 1-513
#$ -tc 513


listName=$(awk -F "\t" "NR==$SGE_TASK_ID {print \$1}" CAD_MAGMA_Global_2_GeneAndSetAnalysis.param)
trait=$(awk -F "\t" "NR==$SGE_TASK_ID {print \$2}" CAD_MAGMA_Global_2_GeneAndSetAnalysis.param)
gwasParam=$(awk -F "\t" "NR==$SGE_TASK_ID {print \$3}" CAD_MAGMA_Global_2_GeneAndSetAnalysis.param)
echo $listName
echo $trait

# create gene set file from each interactor list
mkdir -p temp

tail -n +2 CAD_MasterInteractorTable.txt | \
awk -v var="$listName" 'BEGIN{print var} {if ($1==var && $5=="TRUE") print $4}' | \
tr '\n' '\t' > temp/${listName}.GLOBAL.${trait}.InteractorGeneSet.txt


# directory to store MAGMA output files
mkdir -p magma_output

# MAGMA directory
magmaDir=Software/MAGMA

# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.GLOBAL.genes.annot \
--pval ${gwasParam} \
--out magma_output/${listName}.GLOBAL.${trait}

# Gene Set Analysis Step (calculate enrichment of interactors compared to other genes)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.GLOBAL.${trait}.genes.raw \
--set-annot temp/${listName}.GLOBAL.${trait}.InteractorGeneSet.txt \
--out magma_output/${listName}.GLOBAL.${trait}


# delete temp files
rm temp/${listName}.GLOBAL.${trait}.InteractorGeneSet.txt
