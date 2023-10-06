## libs
setwd('~/Projects/03_MICOM/')
library(dplyr)
library(xlsx)
library(biomaRt)
devtools::load_all('~/Projects/08_genoppi_lagelab/Genoppi/')

# check input
#haarst2017 <- read.csv('genelist/haarst2017/haarst2017_tablexx_gwas_snps.txt', sep = '\t')
#haarst2017 <- haarst2017[complete.cases(haarst2017),]
#haarst2017 <- haarst2017[haarst2017$Gwas.SNP == 1, c(2,3,4,5)]
#colnames(haarst2017) <- c('ID', 'CHR', 'POS', 'P')
#write.table(haarst2017,'genelist/haarst2017/haarst2017_tablexx_gwas_indexsnps.txt', sep = '\t', row.names = F, quote = F)

### extend 160 SNP regions with 50 kb upstream and downstream (anneal generated using plink2 with EUR reference panel, r2 > 0.6)

## merge with P-values
haarst2017 <- read.table('genelist/haarst2017/27AUG20_haarst2017_160snps_tags.tags.list', header = T)
haarst2017_pvalues <- read.csv('genelist/haarst2017/haarst2017_tablexx_gwas_snps.txt', sep = '\t')[,c('Gwas.SNP','Gwas.SNP.1', 'P.value')]
haarst2017_pvalues <- haarst2017_pvalues[complete.cases(haarst2017_pvalues) & haarst2017_pvalues$Gwas.SNP == 1, ]
haarst2017 <- merge(haarst2017, haarst2017_pvalues, by.x = 'SNP', by.y = 'Gwas.SNP.1')

## Prepare data format (colnames, numerics)
haarst2017$region.start <- as.numeric(haarst2017$LEFT) #as.numeric(unlist(lapply(strsplit(haarst2017$left, '\\-'), function(x) x[1])))
haarst2017$region.end <- as.numeric(haarst2017$RIGHT) #as.numeric(unlist(lapply(strsplit(haarst2017$position, '\\-'), function(x) x[2])))
colnames(haarst2017)[2] <- 'chromosome_name'
colnames(haarst2017)[9] <- 'P'
haarst2017$SNP.NOVEL <- as.numeric(haarst2017$SNP %in% read.table('genelist/haarst2017/haarst2017_g1000_eur_anneals.tsv', header = T)$SNP)

# add source of non-novel loci (study name, url etc.)
origin <- read.table('genelist/haarst2017/haarst2017_table2.txt', sep = '\t', header = T)
origin <- origin[,c('SNP','study.name')]
tmp <- data.frame(SNP = unique(origin$SNP))
tmp1 <- lapply(tmp$SNP, function(x) paste(apply(origin[origin$SNP %in% x, ], 1, paste, collapse = ' '), collapse = '; '))
names(tmp1) <- tmp$SNP
origin <- stack(tmp1)[,c(2,1)]
colnames(origin) <- c('SNP','study')
haarst2017 <- merge(haarst2017, origin, by = 'SNP', all.x = T)

# combine data
anneal <- haarst2017
#anneal$TAGS <- NULL

# append 50 kb to each side
rng = 50000
anneal$intervalStart <- anneal$region.start-rng
anneal$intervalEnd <- anneal$region.end+rng

## use biomatr
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","hgnc_symbol","chromosome_name","start_position","end_position", 'gene_biotype'), mart = ensembl)
colnames(mapping) <- c("ensembl_gene_id", "external_gene_name","gene_name", "chromosome_name","geneStart","geneEnd", "gene_type")
ranges <- mapping[mapping$chromosome_name %in% anneal$chromosome_name,]
colnames(ranges)[4:5] <- c('intervalStart','intervalEnd')

## combine by chromosome name, i.e. all combinations of regions each chromosomne, 
## with all possible genes on that chromosome: (snp region x all genes)
ranges <- merge(mapping, anneal, by = 'chromosome_name', all.y = T)

## subet genes
genes <- ranges[with(ranges, 
                     (geneEnd >= intervalStart & geneEnd <= intervalEnd) | # should get partially overlapping genes
                     (geneStart >= intervalStart & geneEnd <= intervalEnd) | # get completely overlapping genes
                     (geneStart <= intervalEnd & geneStart >= intervalStart) | # get partially overlapping genes
                     (intervalStart >= geneStart & intervalEnd <= geneEnd) # get genes that are overlapping the region entirely
),]

# 
#baits <- c('CACNA1C', 'CACNB2', 'CSMD1', 'CUL3', 'GRIN2A', 'HCN1', 'TCF4', 'ZNF804A')
#all(baits %in% genes$gene_name)

# merge with pvalues

# which are single locus genes
genes$single.locus <- genes$SNP %in% names(table(genes$SNP)[table(genes$SNP) == 1])
write.xlsx(genes, file = 'genelist/27AUG20_haarst2017_ensembl_anneal_mapping.xlsx')
write.table(genes, 'genelist/27AUG20_haarst2017_ensembl_anneal_mapping.tsv', row.names = F, quote = F, sep = '\t')

# check with genoppi results
genoppi_genes = unique(unlist(genoppi::get_gene_from_snp(haarst2017$SNP)))
sum(genoppi_genes %in% genes$external_gene_name)/length(genoppi_genes) # ~ 97%

# what baits are present?
baits <- c('ARHGEF26', 'BCAS3', 'JCAD', 'KCNK5', 'FLT1', 'FN1', 'PLPP3', 'KIAA1462', 'HDAC9', 'EDNRA', 'ADAMTS7', 'PPAP2B')
baits %in% genes$gene_name # JCAD ~ KIAA1462. Therefore only PLPP3 is missing.
baits[baits %nin% genes$gene_name]

# write only protein coding gene lists
haarst2017_pc = data.frame(geneName = sort(unique(genes[genes$gene_type == 'protein_coding',]$external_gene_name)))
write.table(haarst2017_pc, 'genelist/27AUG20_haarst2017_protein_coding_genelist.tsv', row.names = F, quote = F, sep = '\t')

