


desc = unlist(lapply(dfs1, function(x) x$Description))
acc = unlist(lapply(dfs1, function(x) x$Accession))

df = data.frame(Accession = acc, Description = desc, stringsAsFactors = F)

df$uniprot.accession.full = accession.uniprot.get(df$Accession)

# get isoforms 
x = uniprot.isoform.get(df$uniprot.accession.full)
df$uniprot.accession.id = x$uniprot
df$uniprot.accession.isoform = x$id

# use table to map to gene id
df$gene.aliases.mapped = map_uniprot_to_gene(as.character(df$uniprot.accession.id))
df$gene.first = unlist(lapply(strsplit(df$gene.aliases.mapped, split = '\\ '), function(x) x[1]))
df$gene.first.approved.symbol = df$gene.first %in% hgnc$Approved.symbol


#df = df[is.na(df$gene.aliases.mapped),]
#df$gene.description = NA
df$gene.description[grepl('GN\\=', df$Description)] = description.keyword(keyword = 'GN=', from = df$Description[grepl('GN\\=', df$Description)])





