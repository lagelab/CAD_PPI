


# 

# dfs1 is a charcater of uniprot IDs
#whitehead.uniprot = unique(unlist(dfs1))


g = unlist(dfs2)
g = unlist(strsplit(g,split = '\\|'))
g = unlist(lapply(strsplit(g,'\\-|\\.'), function(x) x[1] ))
broad.uniprot = unique(g)


x = unique(c(broad.uniprot, whitehead.uniprot))

write.table(x, file = '~/Desktop/uniprot_to_hgnc_whitehead+broad.tsv', quote = F, col.names = F, row.names = F)

