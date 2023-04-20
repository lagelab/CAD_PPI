# expect gnomad LOF table from the official gnomad site

gnomad_dictify <- function(gnomad){
  
  genes = gnomad$gene
  gnomad_genes = lapply(genes, function(gene){
    bool = genes %in% gene
    rows = gnomad[bool, ]
    chr = rows$chromosome
    start_position = rows$start_position
    end_position = rows$start_position
    gene_id = rows$gene_id
    return(data.frame(gene_id = gene_id, chr = chr, start_position = start_position, end_position = end_position))
  })
  names(gnomad_genes) = genes
  return(gnomad_genes)
}