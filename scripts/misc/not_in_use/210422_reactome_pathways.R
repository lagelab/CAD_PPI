# generate a table for reactome pathways

# subset by reactome and setup genesets
data("msigdb_c2_table")
reactome_pathways <- unique(msigdb_c2_table$Set.name)
reactome_pathways <- reactome_pathways [grepl('REACTOME',reactome_pathways )]
msigdb_c2_table$gene <- msigdb_c2_table$Gene.symbol
msigdb_c2_table$significant <- TRUE
reactome <- msigdb_c2_table
reactome <- reactome[reactome$Set.name %in% reactome_pathways,]
write.csv(reactome, 'derived/tables/210422_reactome_msigdb_c2_table.csv', quote = F, row.names = F)
