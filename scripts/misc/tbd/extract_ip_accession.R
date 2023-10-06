

extract_ip_accession <- function(files = list.files('~/Projects/03_MICOM/data/genoppi_input/run030/', full.names = T)){
  
   data = lapply(files, function(f) read.csv(f, sep = '\t'))
   accession_numbers = unlist(lapply(data, function(d) d[,grepl('ccession',colnames(d)) & c(T, F)])) 
   accession_reduced = unique(accession.uniprot.get(accession_numbers))
   accession_split = uniprot.isoform.get(accession_reduced)
   # note some files do not contain acession numbers,
   # since they were given to me this way (gene, rep1, rep2)
   
   return(accession_split$uniprot)
   
}