
# Genoppi workflow for mapping SNPs
# to genes.

load('/Users/flassen/Projects/04_genoppi/Genoppi-master/data/snp_to_gene.RData') # April’s snp_to_gene mapping

snp_to_gene <- function(snp){
  
  if ('genes_snps' %nin% ls(envir = .GlobalEnv)) 
  {write('reading SNP/Gene refrence..',stderr()); data('genes_snps')}
  genes = names(genes_snps)
  result = lapply(genes, function(g){
    table_snps = genes_snps[[g]]
    check = table_snps %in% as.vector(snp) 
    if (any(check)){return(table_snps[check])}
  })
  names(result) = genes
  return(null.omit(result))
}

x= c('rs2820315')
