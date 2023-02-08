##########################################################################################
## UGER job submission script
## run conditional MAGMA analysis for each CAD PPI network (ints vs. non-ints)
## for GWAS traits: CAD, Aorta-AA, Aorta-DA, Stroke-AS, Stroke-AIS, Stroke-LAS,
##                  Stroke-CES, Stroke-SVS, height
##
## Author: Yu-Han Hsu
##########################################################################################

#$ -cwd
#$ -N uger.magma.cond
#$ -l h_vmem=4g
#$ -l h_rt=24:00:00
#$ -t 1-57
#$ -tc 57


listName=$(awk "NR==$SGE_TASK_ID+1 {print \$1}" CAD_Network-Cell-Bait_Mapping.txt)
echo $listName

# create gene.loc and gene set files from each interactor list
mkdir -p temp

tail -n +2 CAD_MasterInteractorTable.txt | \
awk -v var="$listName" '{if ($1==var) {print $4,$6,$7,$8}}' > \
temp/${listName}.gene.loc

tail -n +2 CAD_MasterInteractorTable.txt | \
awk -v var="$listName" 'BEGIN{print var} {if ($1==var && $5=="TRUE") print $4}' | \
tr '\n' '\t' > temp/${listName}.InteractorGeneSet.txt


# directory to store MAGMA output files
mkdir -p magma_output

# MAGMA directory
magmaDir=Software/MAGMA

# Annotation Step (SNP to gene mapping)
${magmaDir}/magma_v1.09_static/magma \
--annotate window=50 \
filter=${magmaDir}/ReferenceData/g1000_eur/g1000_eur.noMHC.bim \
--snp-loc ${magmaDir}/ReferenceData/g1000_eur/g1000_eur.bim \
--gene-loc temp/${listName}.gene.loc \
--out magma_output/${listName}.COND


### CAD
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.COND.genes.annot \
--pval ${magmaDir}/GwasData/CAD/VanDerHarst.CircRes2018.CAD_META.txt snp-id=oldID pval=P-value N=547261 \
--out magma_output/${listName}.COND.CAD

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.COND.CAD.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.COND.CAD


### Aorta-AA
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.COND.genes.annot \
--pval ${magmaDir}/GwasData/Aorta/invnorm_max_aa_diam.maf0.001.tsv snp-id=SNP pval=P_BOLT_LMM N=33420 \
--out magma_output/${listName}.COND.Aorta-AA

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.COND.Aorta-AA.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.COND.Aorta-AA


### Aorta-DA
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.COND.genes.annot \
--pval ${magmaDir}/GwasData/Aorta/invnorm_max_da_diam.maf0.001.tsv snp-id=SNP pval=P_BOLT_LMM N=33420 \
--out magma_output/${listName}.COND.Aorta-DA

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.COND.Aorta-DA.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.COND.Aorta-DA


### Stroke-AS
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.COND.genes.annot \
--pval ${magmaDir}/GwasData/Stroke/MEGASTROKE.1.AS.TRANS.out snp-id=MarkerName pval=P-value N=521612 \
--out magma_output/${listName}.COND.Stroke-AS

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.COND.Stroke-AS.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.COND.Stroke-AS


### Stroke-AIS
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.COND.genes.annot \
--pval ${magmaDir}/GwasData/Stroke/MEGASTROKE.2.AIS.TRANS.out snp-id=MarkerName pval=P-value N=514791 \
--out magma_output/${listName}.COND.Stroke-AIS

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.COND.Stroke-AIS.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.COND.Stroke-AIS


### Stroke-LAS
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.COND.genes.annot \
--pval ${magmaDir}/GwasData/Stroke/MEGASTROKE.3.LAS.TRANS.out snp-id=MarkerName pval=P-value N=461138 \
--out magma_output/${listName}.COND.Stroke-LAS

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.COND.Stroke-LAS.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.COND.Stroke-LAS


### Stroke-CES
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.COND.genes.annot \
--pval ${magmaDir}/GwasData/Stroke/MEGASTROKE.4.CES.TRANS.out snp-id=MarkerName pval=P-value N=463456 \
--out magma_output/${listName}.COND.Stroke-CES

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.COND.Stroke-CES.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.COND.Stroke-CES


### Stroke-SVS
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.COND.genes.annot \
--pval ${magmaDir}/GwasData/Stroke/MEGASTROKE.5.SVS.TRANS.out snp-id=MarkerName pval=P-value N=466160 \
--out magma_output/${listName}.COND.Stroke-SVS

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.COND.Stroke-SVS.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.COND.Stroke-SVS


### height
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.COND.genes.annot \
--pval ${magmaDir}/GwasData/Height/Meta-analysis_Wood_et_al+UKBiobank_2018.txt ncol=N \
--out magma_output/${listName}.COND.height

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.COND.height.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.COND.height


# delete temp files
rm temp/${listName}.gene.loc temp/${listName}.InteractorGeneSet.txt
