
# take an accession number or vector of numbers
# and return their gene name without their flag.tag
get_flag_bait_gene <- function(accession){
  
  res = strsplit(accession, '\\|')
  FLAG = unlist(lapply(res, function(x) {
    bool = grepl('FLAG',x)
    if (any(bool)){
      return(x[bool])
    } else {return(NA)}
  }))
  
  FLAG = gsub('(FLAG)|(3x)|(\\-)|(\\!)|(\\:)', '', FLAG)
  accession[!is.na(FLAG)] <- FLAG[!is.na(FLAG)]
  return(accession)
}


standardize_bait_names <- function(gene){
  
  ## manuel mapping
  # baits
  gene = gsub("^KIAA1462$", 'JCAD', gene)
  
  # non-baits
  gene = gsub("^SSPO$", 'SSPOP', gene) # https://www.genecards.org/cgi-bin/carddisp.pl?gene=SSPO
  gene = gsub("^TWISTNB$", "POLR1F", gene) # https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:18027
  gene = gsub("^HDGFRP3$", "HDGFL3", gene) # https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:24937
  gene = gsub("^CD3EAP$", "POLR1G", gene)
  gene = gsub("^WDR60$", "DYNC2I1", gene)  # https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/21862
  gene = gsub("^ASNA1$", "GET3", gene)
  gene = gsub("^TSTA3$", "GFUS", gene)
  gene = gsub("^LRMP$", "IRAG2", gene)
  gene = gsub("^LOC110117498\\-PIK3R3$", "P3R3URF-PIK3R3", gene)
  gene = gsub("^KIAA0556$", "KATNIP", gene)    

  # read table for automatic mapping
  hgnc.new <- read.csv('data/29MAY20_previous_hgnc_symbol_to_current_hgnc.txt', sep = '\t')
  hgnc.to <- lapply(gene, function(x) unique(as.character(hgnc.new[hgnc.new$Previous.symbol %in% x, 'Approved.symbol']))[1])
  bool = unlist(lapply(hgnc.to, function(x) is.na(x) |identical(x, character(0))))
  gene[!bool] <- unlist(hgnc.to[!bool])
  
  
  
  return(gene)
  
}







