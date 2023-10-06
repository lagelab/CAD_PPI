
if (F){
  ## initial analysis
  .libPaths('~/Toolbox/rlib')
  setwd('~/Projects/03_MICOM/')
  library(pRroteomics)
  library(dplyr)
  library(openxlsx)
  devtools::load_all('~/Toolbox/packages/pRoteomics/')

  write.ip = function(martin, bait, cell, bait.vs = 'Mock', facility = 'Broad', note = 'expanded',imputation = 'MinImputed'){
    return(paste0('Genoppi_',facility,'_martin',martin,'_',imputation,'.',bait,'vs',bait.vs,'.',cell,'.',note,'.tsv'))
  }
}


(files = list.files('~/Projects/03_MICOM/data/raw/martinBroad/', recursive = T, full.names = T, pattern = 'xlsx')) #, pattern = 'nonSGS'))

###############
## ADAMTS7
## cell SMC
# ITRAQ (4plex)
# 114 - Control
# 115 - Control
# 116 - AdamTS7 Lysate IP
# 117 - AdamTS7 Lysate IP
file = "/Users/flassen/Projects/03_MICOM/data/raw/martinBroad//ADAMTS7/20170123_bait_found_HCASMC_only/MacDonald_20170123_HCASMCConditionedMedia+Lysate_DataHandoff_proteinProteinCentricColumnsExport.xlsx"
bait = 'ADAMTS7'
bait_found = TRUE
infile = read.xlsx(file, sheet = 7)
data = infile[,c("groupNum",
                 "subgroupNum",
                 "gene_symbol",
                 "accession_number",
                 "accession_numbers",
                 "entry_name",
                 "LysateIP..Rep1..AdamTS7.vs.IgGControl.(116/114)",
                 "LysateIP..Rep2..AdamTS7.vs.IgGControl.(117/114)",
                 "Gaelen/AdamTS7/Jan2017/Lysate_Longrun.unique_peptides")]
colnames(data) = c("groupNum","subgroupNum",'id', 'accession_number', 'accession_numbers', 'description', 'rep1', 'rep2', 'unique.peptides')

#data[data$unique.peptides == 1,] # only the bait is found once, no need for further filtering
#expand_ids = expand_accession_id(data$accession)
#data$gene = uniprot_to_hgnc(expand_ids$uniprot)
#data %>% normalize() %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotScatter(bait, title = basename(file))
#data %>% normalize() %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotVolcano(bait, title = basename(file))
#data %>% normalize() %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotOverlap(bait, reference = interactors('ADAMTS7'),  title = basename(file))
write.table(data, paste0('data/genoppi_input/',write.ip('Broad', bait = bait, cell='SMC')), quote = F, sep = '\t')

#int = interactors('ADAMTS7')
#int[int$significant, ]$gene %in% data$gene


# Note: No InWeb interactors could be found. Futhermore, the bait
# is only associated with one unique peptide. Variance is also not sufficient
# to calculate moderated ttest.
################

###############
## ADAMTS7
## cell EC
# ITRAQ (4plex)
# 114 - Control
# 115 - Control
# 116 - AdamTS7 Lysate IP
# 117 - AdamTS7 Lysate IP
file = "/Users/flassen/Projects/03_MICOM/data/raw/martinBroad//ADAMTS7/20170815_bait_found_HAEC_only/Zhu_20170815_HAEClysate_HCASMClysate_Combined_modtdev.xlsx"
stopifnot(file.exists(file))
bait = 'ADAMTS7'
bait_found = TRUE
infile = read.xlsx(file, sheet = 5)
data = infile[,c("subgroupNum",
                 "id",
                 "id.mapped",
                 "accession_numbers",
                 'entry_name',
                 "EC_LongRun..114.116..ADAMTS7_r1.EV_r1",
                 "EC_LongRun..115.117..ADAMTS7_r2.EV_r2",
                 "EC_LongRun.unique_peptides",
                 'species')]
colnames(data) = c("subgroupNum",'accession_number','id', 'accession_numbers', 'description', 'rep1', 'rep2', 'unique.peptides','species')
#expand_ids = expand_accession_id(data$id)
#expand_ids$hgnc = uniprot_to_hgnc(expand_ids$uniprot)
#data$id = expand_ids # verify mapping.
#data$gene[is.na(data$gene)] = as.character(data$id$uniprot[is.na(data$gene)])
#nrow(data[data$unique.peptides < 2,]) # 274 non-mapped peptides
#data = data[data$unique > 1 ,c('gene', 'accession', 'rep1', 'rep2')]
#data = impute.gaussian(data); sum(data$imputed)
#data %>% normalize() %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotScatter('ADAMTS7') # 0.2 correlation
#data %>% normalize() %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotVolcano('ADAMTS7') # nothing is significant
write.table(data, paste0('data/genoppi_input/',write.ip('Broad', bait = bait, cell='EC')), quote = F, sep = '\t')

# note: data seems non-significant. Low scatter plot correlations between replicates.
# Probably failed IP.
########


###############
## EDN1
## cell EC and SMC combined
# ITRAQ (4plex)
# 114 - Control EC SMC
# 115 - Control EC SMC
# 116 - EDN1 IP EC SMC
# 117 - EDN1 IP EC SMC
file = "/Users/flassen/Projects/03_MICOM/data/raw/martinBroad//EDN1/20170221_bait_found/20170221_EDN1_SMC_and_EC_Combined_nSGS.xlsx"
bait = 'EDN1'
bait_found = TRUE
infile_EC = read.xlsx(file, sheet = 5)
infile_SMC = read.xlsx(file, sheet = 6)

# EC
data_EC = infile_EC[,c("Gene_Symbol",
                 "accession_number",
                 "accession_numbers",
                 'entry_name',
                 "EC.(116/114).EDN1_R1/Control_R1",
                 "EC.(117/115).EDN1_R2/Control_R2",
                 "EC.unique_peptides")]
colnames(data_EC) = c('id', 'accession_number', 'accession_numbers', 'description', 'rep1', 'rep2', 'unique.peptides')
#data_EC_expanded = expand_accession_id(data_EC$accession)
#data_EC$gene = uniprot_to_hgnc(data_EC_expanded$uniprot)
#data_EC = data_EC[data_EC$unique.peptides > 1, c('gene', 'accession', 'rep1', 'rep2')]
#data_EC = impute.gaussian(data_EC)
#data_EC %>% normalize() %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotScatter('EDN1', title = "Broad (EDN1) SMC")
#data_EC %>% normalize() %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotVolcano('EDN1', title = "Broad (EDN1) EC")
write.table(data_EC, paste0('data/genoppi_input/',write.ip('Broad', bait = bait, cell='EC')), quote = F, sep = '\t')

# SMC
data_SMC = infile_SMC[,c("Gene_Symbol",
                         "accession_number",
                         "accession_numbers",
                         'entry_name',
                       "SMC.(114/116).EDN1_R1/Control_R1",
                       "SMC.(115/117).EDN1_R2/Control_R2",
                       "SMC.unique_peptides")]
colnames(data_SMC) = c('id', 'accession_number', 'accession_numbers', 'description', 'rep1', 'rep2', 'unique.peptides')
#bait %in% data_EC$gene; bait %in% data_SMC$gene # bait is there
#data_SMC_expanded = expand_accession_id(data_SMC$accession)
#data_SMC$gene = uniprot_to_hgnc(data_SMC_expanded$uniprot)
#data_SMC = data_SMC[data_SMC$unique.peptides > 1, c('gene', 'accession', 'rep1', 'rep2')]
#data_SMC = impute.gaussian(data_SMC)
#data_SMC %>% normalize() %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotScatter('EDN1', title = "Broad (EDN1) SMC")
#data_SMC %>% normalize() %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotVolcano('EDN1', title = "Broad (EDN1) EC")
write.table(data_SMC, paste0('data/genoppi_input/',write.ip('Broad', bait = bait, cell='SMC')), quote = F, sep = '\t')

# note: EDN1 in SMC is an OK IP.
# whereas the EC IP is not useable.
########

###############
## HDAC9
## cell EC and SMC combined
# ITRAQ (4plex)
# 114 - MOCK
# 115 - MOCK
# 116 - Flag-HDAC9
# 117 - Flag-HDAC9
file = "/Users/flassen/Projects/03_MICOM/data/raw/martinBroad//HDAC9/20170511_HDAC9_FlagSearched_DataHandoffproteinProteinCentricColumnsExport.1.xlsx"
bait = 'HDAC9'
bait_found = TRUE
infile = read.xlsx(file, sheet = 7)
data = infile[,c("subgroupNum",
                 "geneName",
                 "accession_number",
                 "accession_numbers",
                 'entry_name',
                 "Gaelen/HDAC9/HDAC9_HCASMC_BRP_FlagSearch.unique_peptides",
                 "HCASMC..MedianNormalized.Log2..Rep1..Flag.HDAC9.vs.Mock..116.v.114",
                 "HCASMC..MedianNormalized.Log2..Rep2..Flag.HDAC9.vs.Mock..117.v.115",
                 "Gaelen/HDAC9/HDAC9_HAEC_BRP_FlagSearch.unique_peptides",
                 "HAEC..MedianNormalized.Log2..Rep1..Flag.HDAC9.vs.Mock..116.v.114",
                 "HAEC..MedianNormalized.Log2..Rep2..Flag.HDAC9.vs.Mock..117.v.115")]
colnames(data) = c('subgroupNum', 'id','accession_number', 'accession_numbers', 'description', 'smc.unique.peptides', 'smc.rep1', 'smc.rep2', 'ec.unique.peptides', 'ec.rep1', 'ec.rep2')
#data_expanded = expand_accession_id(data$accession)
#data$gene = as.character(uniprot_to_hgnc(data_expanded$uniprot))


data_smc = data[, c("subgroupNum", 'id', 'accession_number','accession_numbers', 'description', 'smc.rep1', 'smc.rep2','smc.unique.peptides')]
data_ec = data[, c("subgroupNum", 'id', 'accession_number','accession_numbers', 'description', 'ec.rep1', 'ec.rep2','ec.unique.peptides')]
#colnames(data_smc) = c('gene','accession', 'rep1', 'rep2')
#colnames(data_ec) = c('gene','accession', 'rep1', 'rep2')
#data_ec = impute.gaussian(data_ec)
#data_ec %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotScatter(bait)
#data_ec %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotVolcano(bait)
write.table(data_ec, paste0('data/genoppi_input/',write.ip('Broad', bait = bait, cell='EC')), quote = F, sep = '\t')


#data_smc = impute.gaussian(data_smc)
#data_smc %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotScatter(bait)
#data_smc %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotVolcano(bait)
write.table(data_smc, paste0('data/genoppi_input/',write.ip('Broad', bait = bait, cell='SMC')), quote = F, sep = '\t')

# note: both datasets have FDR > 0.1, even though two isoforms
# of the bait is enriched, the data qaulity is very low.
######

###############
## PHACTR1
## cell EC and SMC combined (bait only found in EC)
# ITRAQ (4plex)

file = "/Users/flassen/Projects/03_MICOM/data/raw/martinBroad//PHACTR1/20170203_bait_found_EC_only/20170131_PHACTR1_EC_and_SMC_nSGS.xlsx"
bait = 'PHACTR1'
bait_found = TRUE
infile = read.xlsx(file, sheet = 7)
data = infile[,c("subgroupNum",
                 "ID",
                 "accession_number",
                 "accession_numbers",
                 "entry_name",
                 "HCASMC.(114/116).PHACTR1_R1/Cntl_R1",
                 "HCASMC.(115/117).PHACTR1_R2/Cntl_R2",
                 "HAEC.(116/114).PHACTR1_R1/Cntl_R1",
                 "HAEC.(117/115).PHACTR1_R2/Cntl_R2",
                 "HAEC.unique_peptides",
                 "HCASMC.unique_peptides")]
colnames(data) = c("subgroupNum",'id','accession_number', 'accession_numbers', 'description', 'smc.rep1', 'smc.rep2', 'ec.rep1', 'ec.rep2', 'ec.unique.peptides', 'smc.unique.peptides')

#data_expanded = expand_accession_id(data$accession)
#data$gene = as.character(uniprot_to_hgnc(data_expanded$uniprot))
data_smc = data[, c("subgroupNum",'id','accession_number', 'accession_numbers', 'description', 'smc.rep1', 'smc.rep2','smc.unique.peptides')]
data_ec = data[, c("subgroupNum",'id','accession_number', 'accession_numbers', 'description','ec.rep1', 'ec.rep2','ec.unique.peptides')]
colnames(data_smc) = c("subgroupNum",'id','accession_number', 'accession_numbers', 'description', 'rep1', 'rep2','unique.peptides')
colnames(data_ec) = c("subgroupNum",'id','accession_number', 'accession_numbers', 'description', 'rep1', 'rep2','unique.peptides')
#data_ec = impute.gaussian(data_ec)
#data_ec %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotScatter(bait)
#data_ec %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotVolcano(bait)
write.table(data_ec, paste0('data/genoppi_input/',write.ip('Broad', bait = bait, cell='EC')), quote = F, sep = '\t')


#data_smc = impute.gaussian(data_smc)
#data_smc %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotScatter(bait)
#data_smc %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotVolcano(bait)
write.table(data_smc, paste0('data/genoppi_input/',write.ip('Broad', bait = bait, cell='SMC')), quote = F, sep = '\t')


# note: EC bait is found, but very low correlations, wheras for SMC
# bait is not found but good correlations.
#########


###############
## PLPP3
## cell EC and SMC combined (bait only found in EC)
# ITRAQ (4plex)
file = "/Users/flassen/Projects/03_MICOM/data/raw/martinBroad//PLPP3/20170306_bait_found/20170306_PLPP3_SMC_and_EC_nSGS.xlsx"
bait = 'PLPP3' #(UNIPROT_ID=O14495)
bait_found = TRUE
infile = read.xlsx(file, sheet = 8)
data = infile[,c("subgroupNum",
                 "Gene_Symbol",
                 "accession_number",
                 "accession_numbers",
                 "EC.(114/116).PLPP3_R1/Control_R1",
                 "EC.(115/117).PLPP3_R2/Control_R2",
                 "SMC.(116/114).PLPP3_R1/Control_R1",
                 "SMC.(117/115).PLPP3_R2/Control_R2",
                 "EC.unique_peptides","SMC.unique_peptides")]
colnames(data) = c("subgroupNum", "Gene_Symbol", "accession_number","accession_numbers", 'ec.rep1', 'ec.rep2', 'smc.rep1', 'smc.rep2', 'ec.unique.peptides', 'smc.unique.peptides')
#data_expanded = expand_accession_id(data$accession)
#data$gene = as.character(uniprot_to_hgnc(data_expanded$uniprot))
data_smc = data[, c("subgroupNum", "Gene_Symbol", "accession_number","accession_numbers", 'smc.rep1', 'smc.rep2','smc.unique.peptides')]
data_ec = data[, c("subgroupNum", "Gene_Symbol", "accession_number","accession_numbers", 'ec.rep1', 'ec.rep2','ec.unique.peptides')]
colnames(data_smc) = c("subgroupNum", "id", "accession_number","accession_numbers", 'rep1', 'rep2', 'unique.peptides')
colnames(data_ec) = c("subgroupNum", "id", "accession_number","accession_numbers", 'rep1', 'rep2', 'unique.peptides')
#data_ec = impute.gaussian(data_ec)
#data_ec %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotScatter(bait)
#data_ec %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotVolcano(bait)
write.table(data_ec, paste0('data/genoppi_input/',write.ip('Broad', bait = bait, cell='EC')), quote = F, sep = '\t')


#data_smc = impute.gaussian(data_smc)
#data_smc %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotScatter(bait)
#data_smc %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotVolcano(bait)
write.table(data_smc, paste0('data/genoppi_input/',write.ip('Broad', bait = bait, cell='SMC')), quote = F, sep = '\t')


###############
## ARHGEF26 ##

(files = list.files('~/Projects/03_MICOM/data/raw/martinBroad/ARHGEF26/', recursive = T, full.names = T, pattern = 'SGS')) #, pattern = 'nonSGS'))
for (fname in files){

  # import data
  print(fname)
  data_import = read.csv(fname, sep = '\t', stringsAsFactors = F)
  data = data_import
  data$rep1 <- as.numeric(data$rep1)
  data$rep2 <- as.numeric(data$rep2)

  # deal with weird NAs
  index = unique(c(which(is.na(data$rep1)), which(is.na(data$rep2))))
  print(data_import[index, ])
  data <- data[complete.cases(data), ]


  # get cell and bait
  type = gsub('\\.csv','',unlist(strsplit(basename(fname),'_'))[3])
  cell = gsub('','',unlist(strsplit(basename(fname),'_'))[2])
  bait = unlist(strsplit(basename(fname),'_'))[1]

  # logfold change
  print(nrow(data))
  data = data[,c('gene', 'rep1', 'rep2')]
  write.table(data, paste0('data/genoppi_input/',write.ip('Broad', bait=bait, cell=cell, note = type)), quote = F, sep = '\t')


  #calculated = data %>% mttest()
  #print(fname)
  #calculated %>% designate(FDR < 0.1, logFC > 0) %>% plotVolcano('ARHGEF26', title = fname)
  #Sys.sleep(3)

}



## OLD

###############
## section 1 ##
#file = "/Users/flassen/Projects/03_MICOM/data/raw/martinBroad//ADAMTS7/20170123_bait_found_HCASMC_only/MacDonald_20170123_HCASMC_ConditionedMedia+Lysate_DataHandoff_proteinProteinCentricColumnsExport.2.xlsx"
#bait = 'PLPP3'
#bait_found = TRUE
#data = read.xlsx(file, sheet = 5)

# get relevant columns
#cols_itraq = grepl('iTRAQ_[0-9]{3}_total', colnames(data))
#cols_itraq_lysate = grepl('Lysate.+iTRAQ_[0-9]{3}_total', colnames(data))
#cols_itraq_total = cols_itraq & !cols_itraq_lysate & !grepl('ConditionedMedia', colnames(data))
#cols_itraq = grepl('iTRAQ_[0-9]{3}_total', colnames(data))
#cols_unique_peptides = grepl('Lysate_Longrun.unique_peptides', colnames(data))
#cols_itraq_conditioned_media = cols_itraq & !cols_itraq_lysate & grepl('ConditionedMedia', colnames(data))
#cols_acession = grepl('accession_number$', colnames(data))
#data_lysate = data[,cols_acession | cols_unique_peptides | cols_itraq_lysate]
#data_total = data[,cols_acession | cols_unique_peptides | cols_itraq_total]

# conditioned media
#columns = c("Gaelen/AdamTS7/Jan2017/ConditionedMedia_50percInject.iTRAQ_114_total",
#            "Gaelen/AdamTS7/Jan2017/ConditionedMedia_50percInject.iTRAQ_116_total",
#
#)

# Lysate
#columns = c('accession_number', 'Gaelen/AdamTS7/Jan2017/Lysate_Longrun.iTRAQ_114_total',
#            "Gaelen/AdamTS7/Jan2017/Lysate_Longrun.iTRAQ_116_total",
#            "Gaelen/AdamTS7/Jan2017/Lysate_Longrun.iTRAQ_115_total",
#            "Gaelen/AdamTS7/Jan2017/Lysate_Longrun.iTRAQ_117_total")
#data = prepare('PLPP3', data_lysate, cols = columns, filter = NULL, raw=T)
#data$data[, c('gene', 'rep1', 'rep2')] %>% mttest() %>% designate(FDR < 0.1, logFC > 0) %>% plotVolcano(bait)

# Total
#columns = c('accession_number', 'iTRAQ_114_total', 'iTRAQ_116_total', 'iTRAQ_115_total', 'iTRAQ_117_total')
#data = prepare('PLPP3', data_total, cols = columns, filter = NULL, raw=T)
#data$data[, c('gene', 'rep1', 'rep2')] %>% mttest() %>% designate(FDR < 0.1, logFC) %>% plotVolcano(bait)


###############
## section 2 ##
#file = "/Users/flassen/Projects/03_MICOM/data/raw/martinBroad//ADAMTS7/20170123_bait_found_HCASMC_only/MacDonald_20170123_HCASMCConditionedMedia+Lysate_DataHandoff_proteinProteinCentricColumnsExport.xlsx"
#bait = 'PLPP3'
#bait_found = TRUE
#data = read.xlsx(file, sheet = 5)
#colnames(data)
