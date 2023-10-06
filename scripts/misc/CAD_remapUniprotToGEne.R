## remap gene ids based on uniprot website https://www.uniprot.org/mapping/M202002116746803381A1F0E0DB47453E0216320D77A676G.tab

# Load data
files = list.files('data/genoppi_input/run001/', full.names = T)
data = lapply(files, function(df) {read.csv(df,sep='\t')})
names(data) = files


# get list of accession numbers to generate gene mappings
#acc = unique(unlist(lapply(data, function(x) strsplit(as.character(x$accession),'\\|'))))
#acc = expand_accession_id(acc)
#write.csv(data.frame(uniprot=acc$uniprot), '11FEB2020_micom_found_uniprot_ids.csv', row.names = F, quote = F)
mappings = read.table('11FEB2020_uniprot_to_gene_table.txt', header = T)
mappings$From <- as.character(mappings$From)
mappings$To <- as.character(mappings$To)

# function for mapping the uniprot accession numbers
# to the correct gene id specified by uniprot-table
micom_uniprot_to_gene <- function(vec){
  lapply(as.character(vec), function(x){ifelse(x %in% mappings$From, mappings[mappings$From == x,]$To, NA)})
}


if (1==1){
new_data = lapply(data, function(x){

  if ('accession' %in% colnames(x)){
    accession = unlist(lapply(strsplit(as.character(x$accession),'\\|'), function(y) y[1]))
    accession_expanded = expand_accession_id(accession)
    accession_uniprot = as.character(accession_expanded$uniprot)
    print(paste(length(accession_uniprot), '-', length(x$gene)))
    x$ngene = unlist(micom_uniprot_to_gene(accession_uniprot))
  } else {
    x$ngene = x$gene
  }
  return(x)
})
}

names(data) = basename(files)

count = 1
new_data_unique = lapply(new_data, function(x){

#baits = as.character(unique(read.csv('10FEB2020_masterlist_analysis_plan.csv')$Bait))
#baits = baits[baits != '']

#for (x in new_data){

  #browser()
  #if (any(is.na(x$gene))) browser()
  print(names(new_data)[count])
  
  # Map missing ids with 
  df = x %>% mttest(keep=c('imputed','accession','ngene')) %>% designate(FDR < 0.1, logFC > 0)
  df$gene <- as.character(df$gene)
  df$ngene <- as.character(df$ngene)
  #print(head(df))
  #Sys.sleep(60)
  
  
  #print(sum(na.omit(df$gene == df$ngene))/nrow(df))
  
  # re-map missing genes
  print(paste('trying re-mapping', sum(is.na(df$gene)),'NAs (genes) using uniprot-table.' ))
  df[is.na(df$gene), ]$gene = as.character(df[is.na(df$gene), ]$ngene)
  
  if (any(grepl('FLAG', toupper(df$gene)))) print(df[grepl('FLAG', toupper(df$gene)), ])
  
  # re-map names with flag tag
  if (any(grepl('^3xFLAG-ARHGEF26$', df$gene))) df[df$gene == '3xFLAG-ARHGEF26',]$gene <- 'ARHGEF26'
  if (any(grepl('^3xFLAG-ARHGEF26-MT$', df$gene))) df[df$gene == '3xFLAG-ARHGEF26-MT',]$gene <- 'ARHGEF26'
  if (any(grepl('^3xFLAG-ARHGEF26-WT$', df$gene))) df[df$gene == '3xFLAG-ARHGEF26-WT',]$gene <- 'ARHGEF26'
  if (any(grepl('JCAD-3xFLAG', df$gene)))  df[df$gene == 'JCAD-3xFLAG',]$gene <- 'JCAD'
  if (any(grepl('KIAA1462', toupper(df$gene)))) df[df$gene == 'KIAA1462',]$gene <- 'JCAD'
  if (any(grepl('EDN1-3xFLAG', df$gene)))  df[df$gene == 'EDN1-3xFLAG',]$gene <- 'EDN1'
  if (any(grepl('PHACTR1-FLAG', df$gene)))  df[df$gene == 'PHACTR1-FLAG',]$gene <- 'PHACTR1'
  if (any(grepl('!FINC', df$gene)))  df[df$gene == '!FINC',]$gene <- 'FN1'
  if (any(grepl('3xFLAG-BCAS3', df$gene)))  df[df$gene == '3xFLAG-BCAS3',]$gene <- 'BCAS3'
  if (any(grepl('KCNK5-3xFLAG', df$gene)))  df[df$gene == 'KCNK5-3xFLAG',]$gene <- 'KCNK5' 
  if (any(grepl('ADAMTS7-3xFLAG', df$gene)))  df[df$gene == 'ADAMTS7-3xFLAG',]$gene <- 'ADAMTS7'
  if (any(grepl('3xFLAG-HDAC9', df$gene)))  df[df$gene == '3xFLAG-HDAC9',]$gene <- 'HDAC9'
  
  
  # need to handle FLT1 - JCAD IP individually
  
  if (any(df$ngene %in% baits)) print(df[df$ngene %in% baits,])
  
  # handle special cases
  if (any(grepl('FLAG', toupper(df$gene)))) browser()
  #if (any(grepl('WT', toupper(df$gene)) & !grepl('WTAP', toupper(df$gene)) )) browser()
  #if (any(grepl('MT', toupper(df$gene)))) browser()
  
  ## need to decide on a mapping
  
  #print(df[is.na(df$ngene), ])
  
  # Remove duplicates
  #if (any(duplicated(df$gene))) browser()
  tab = sort(table(df$gene),decreasing = T)
  # find non-unique rows (based on gene name) and only keep
  # ones that have the lowest FDR
  rm = lapply(rownames(tab[tab > 1]), function(gene_name){
    subdf = df[df$gene %in% gene_name, ]
    subdf = subdf[order(subdf$FDR), ]
    #print(subdf)
    rmdf = subdf[2:nrow(subdf), ]
    #print(rownames(rmdf))
    return(rownames(rmdf))
  })
  
  # remove filtered rownames
  rm = unlist(rm)
  df = df[rownames(df) %nin% rm, ]
  
  count = count + 1
  #print('OK!')
  return(df)
})


# write data
lapply(names(new_data_unique), function(file){
  outfile = gsub('run001','run004', file)
  df = new_data_unique[[file]]
  df$ngene = NULL
  df$accession = NULL
  df$significant = NULL
  write.table(df, outfile, quote = F, sep = '\t')
})





#missingness_biomart = round(100*(sum(is.na(x$gene))/length(x$gene)),2)
#missingness_uniprot = round(100*(sum(is.na(x$ngene))/length(x$gene)),2)
#print(paste('biomart = ',missingness_biomart,' and uniprot =', missingness_uniprot))
#return(ifelse(missingness_uniprot > missingness_biomart, 'biomart','uniprot'))


