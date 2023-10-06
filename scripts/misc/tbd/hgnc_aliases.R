

# get the gene and all its aliases
hgnc_aliases <- function(hgnc){
  
  #hgnc <- read.delim(url("https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_name_aliases&col=gd_pub_ensembl_id&col=md_prot_id&col=gd_pub_eg_id&col=gd_prev_name&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"))
  
  tabl = list()
  for (i in 1:nrow(hgnc)){
    
    # get approved symbols and aliases
    approved.symbol = as.character(hgnc$Approved.symbol[i])
    aliases = unique(unlist(lapply(hgnc$Alias.symbols[i], function(x) strsplit(as.character(x), split = ', '))))
    tmp_tabl = list()
    alias = approved.symbol
    
    if (!identical(aliases, character(0))){
      # each alias will point towards the approved gene name
      for (alias in aliases){
        tmp_tabl[[alias]] = c(approved.symbol, alias)
      }
    }
    
    tmp_tabl[[alias]] = c(approved.symbol, approved.symbol)
    
    # combine minor tables
    df = as.data.frame(do.call(rbind, tmp_tabl))
    colnames(df) = c('from', 'to')
    tabl[[i]] = df
    
  }
  
  # get list of tables
  return(tabl)
  
}



hgnc_symbols <- function(hgnc){
  
  #hgnc <- read.delim(url("https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_name_aliases&col=gd_pub_ensembl_id&col=md_prot_id&col=gd_pub_eg_id&col=gd_prev_name&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"))
  
  tabl = list()
  for (i in 1:nrow(hgnc)){
    
    # get approved symbols and aliases
    approved.symbol = as.character(hgnc$Approved.symbol[i])
    aliases = unique(unlist(lapply(hgnc$Previous.symbols[i], function(x) strsplit(as.character(x), split = ', '))))
    tmp_tabl = list()
    alias = approved.symbol
    
    if (!identical(aliases, character(0))){
      # each alias will point towards the approved gene name
      for (alias in aliases){
        tmp_tabl[[alias]] = c(approved.symbol, alias)
      }
    }
    
    tmp_tabl[[alias]] = c(approved.symbol, approved.symbol)
    
    # combine minor tables
    df = as.data.frame(do.call(rbind, tmp_tabl))
    colnames(df) = c('from', 'to')
    tabl[[i]] = df
    
  }
  
  # get list of tables
  return(tabl)
  
}




