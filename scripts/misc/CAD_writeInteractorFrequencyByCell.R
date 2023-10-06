setwd('~/Projects/03_MICOM/')
devtools::load_all('~/Projects/08_genoppi_lagelab/Genoppi/')
library(xlsx)

# get data (downloaded on 14 AUG 2020 From google drive: 200813_MICOM_InteractorFrequency.xlsx - EC)
data = list(
  EC = read.csv('~/Desktop/200813_MICOM_InteractorFrequency.xlsx - EC.tsv', sep = '\t'),
  SMC = read.csv('~/Desktop/200813_MICOM_InteractorFrequency.xlsx - SMC.tsv', sep = '\t'),
  ECSMC = read.csv('~/Desktop/200813_MICOM_InteractorFrequency.xlsx - EC-SMC.tsv', sep = '\t')
)

# get all baits
all_baits = unique(unlist(lapply(data, function(x) unlist(strsplit(x$Baits, split = ', ')))))

# (A) Get different data sources

# (1) Get all, high-confidence and gold-standard interactors
inweb_dict <- lapply(all_baits, function(bait){
  if (bait == 'KIAA1462') bait <- 'KIAA1462'
  A = get_inweb_list(bait) 
  A = A$gene[A$significant]
  H = get_inweb_list(bait, 'hc')
  H = H$gene[H$significant]
  G = get_inweb_list(bait, 'gs') 
  G = G$gene[G$significant]
  return(list(all = A, hc = H, gs = G))
  })
names(inweb_dict) <- all_baits


# (2_ get roselli genes
roselli <- as.vector(unlist(read.csv('data/micom_roselli.txt', header = F)))

## get haarst GWAS
haarst2017 <- read.table('genelist/27AUG20_haarst2017_ensembl_anneal_mapping.tsv', sep = '\t',stringsAsFactors = F, header = T)
haarst2017 <- haarst2017[haarst2017$gene_type == 'protein_coding',]
c('KIAA1462', 'JCAD') %in% haarst2017$external_gene_name
c('PLPP3', 'PPAP2B') %in% haarst2017$external_gene_name
haarst2017$external_gene_name[haarst2017$external_gene_name == 'KIAA1462'] <- 'JCAD'
haarst2017$external_gene_name[haarst2017$external_gene_name == 'PPAP2B'] <- 'PLPP3'
haarst2017 <- haarst2017[,c('external_gene_name','SNP', 'P.value', 'study')]
haarst2017$study[is.na(haarst2017$study)] <- 'Haarst 2017 (Novel loci)'
colnames(haarst2017) <- c('Gene', 'SNP', 'P', 'study')


# (3) Nelson GWAS

## deal with exome chip variants
nelson.gwas <- read.csv('genelist/Nelsen2017_CAD_GWAS_Supplementary_table_4.csv')
nelson.gwas$exomechip.variants <- grepl('\\*',nelson.gwas$Markername)
nelson.gwas$Markername <- gsub('\\*', '', nelson.gwas$Markername)

## subset by genome wide significant variants
nelson.gwas <- nelson.gwas[nelson.gwas$Pvalue <= 5e-08 & complete.cases(nelson.gwas$Pvalue), ]

## Use Genoppi to do SNP to gene mapping
nelson.snp <- stack(lapply(get_gene_from_snp(nelson.gwas$Markername), function(x) paste(x, collapse = ' ')))
colnames(nelson.snp) <- c('gene','Markername')
nelson.snp <- nelson.snp[complete.cases(nelson.snp$gene),]

## merge with original table
nelson.gwas.mapped <- merge(nelson.gwas, nelson.snp)
nelson.gwas.genes <- unique(unlist(lapply(nelson.gwas.mapped$gene, function(x) unlist(strsplit(x, '\\ ')))))


# (B) combine data 
newdata <- lapply(data, function(df){
  
  # pythonic way of iterating through the data
  tmpdata = lapply(1:nrow(df), function(i){
  #for (i in 1:nrow(df)){
    
    ## get relevant data
    interactor <- df$Interactor[i]
    the_baits <- unlist(strsplit(df$Bait[i], split = '\\, '))
    interactor_plus_baits <- c(interactor, the_baits)
    
    ## check for inweb
    inweb_all = unlist(lapply(the_baits, function(x) interactor %in% inweb_dict[[x]]$all))
    inweb_hc = unlist(lapply(the_baits, function(x) interactor %in% inweb_dict[[x]]$hc))
    inweb_gs = unlist(lapply(the_baits, function(x) interactor %in% inweb_dict[[x]]$gs))
    inweb_df = data.frame(all = inweb_all, hc = inweb_hc, gs = inweb_gs)
    rownames(inweb_df) = the_baits
    
    if (any(unlist(inweb_df))){
      inweb_result = apply(inweb_df,1,sum)
      inweb_result = inweb_result[inweb_result != 0]
      inweb_string = gsub(1,'ALL',gsub(2,'HC',gsub(3, 'GS', inweb_result)))
      inweb_result = paste0(paste0(names(inweb_result),':',interactor , ' (',inweb_string,')', sep = ''), collapse = ' ')
    } else {
      inweb_result = ''
    }
    
    # Check for Roselli
    roselli_in_data = interactor %in% roselli
    if (any(roselli_in_data)){
      roselli_result = paste(interactor[interactor %in% roselli], collapse = ' ')
    } else {
      roselli_result = ''
    }

    # Check for haarst
    haarst_in_data = interactor %in% haarst2017$Gene
    if (any(haarst_in_data)){
      haarst_gene = paste(interactor[interactor %in% haarst2017$Gene], collapse = ' ')
      haarst_all = haarst2017[haarst2017$Gene %in% interactor,]
      haarst_snp = paste(haarst_all$SNP, collapse = ' ')
      haarst_pvalue = paste(haarst_all$P, collapse = ' ')
      haarst_study = paste(haarst_all$study, collapse = ' ')
      
    } else {
      haarst_gene = ''
      haarst_snp = ''
      haarst_pvalue = ''
      haarst_study = ''
    }
    
    # Check for Nelson GWAS
    nelson_in_data = interactor %in% nelson.gwas.genes
    if (any(nelson_in_data)){
      nelson_genes = paste(interactor[nelson_in_data], collapse = ' ')
      nelson_row = unlist(lapply(interactor, function(x){
        splitted = strsplit(nelson.gwas.mapped$gene, split = '\\ ')
        myindex = NA
        for (i in 1:length(splitted))if (x %in% splitted[[i]]) myindex = i
        return(myindex)
      }))
      nelson_row = na.omit(nelson_row)
      selected = do.call(rbind, lapply(nelson_row, function(row) nelson.gwas.mapped[row, ]))
      
      nelson_result_variant = paste(selected$Markername, collapse = ' ')
      nelson_result_pvalue = paste(selected$Pvalue, collapse = ' ')
      nelson_result_annotation = paste(selected$annotation, collapse = ' ')
      nelson_result_qvalue = paste(selected$FDR.Qvalue, collapse = ' ')
      nelson_result_OR1 = paste(selected$OR..95..CI., collapse = ' ')
    } else {
      nelson_genes = ''
      nelson_result_variant = ''
      nelson_result_pvalue = ''
      nelson_result_annotation = ''
      nelson_result_qvalue = ''
      nelson_result_OR1 = ''
    }
    
    # combine result in row
    #df[i, ]
    newrow = data.frame(inweb = inweb_result, roselli = roselli_result, 
                        haarst2017_gene = haarst_gene, haarst2017_snp = haarst_snp, 
                        haarst2017_pvalue = haarst_pvalue, haarst2017_study = haarst_study,
                        nelson.genes = nelson_genes, nelson.variant = nelson_result_variant,
                        nelson.pvalue = nelson_result_pvalue, nelson_qvalue = nelson_result_qvalue,
                        nelson_reult_annotation = nelson_result_annotation, nelson_result_OR = nelson_result_OR1)
    
    combined = cbind(df[i, ], newrow)
    
  })
  
  return(as.data.frame(do.call(rbind, tmpdata)))

})


outstatfile = '27AUG20_MICOM_InteractorFrequencyStatistics.EC.SMC.xlsx'
write.xlsx(newdata$EC, file=outstatfile, sheetName="EC", append = FALSE, row.names=FALSE)
write.xlsx(newdata$SMC, file=outstatfile, sheetName="SMC", append = TRUE, row.names=FALSE)
write.xlsx(newdata$ECSMC, file=outstatfile, sheetName="EC-SMC", append = TRUE, row.names=FALSE)






