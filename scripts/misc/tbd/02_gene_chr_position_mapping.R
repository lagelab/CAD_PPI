setwd('~/Projects/03_MICOM/')


# refined replicates 
paths = list.files('data/genoppi_input/run031', full.names = T) # old: 007
paths.whitehead = !grepl('Broad', paths)
paths.broad = grepl('Broad', paths)

dfs = lapply(paths, function(path){
  df = read.csv(path, sep = '\t')
  return(as.character(df$gene.final))
})

#------------------------------------
# what genes need to be mapped

genes = unique(unlist(dfs))

#-----------------------------
# map the genes to regions

genomic.mapping <- read.csv('data/02JUN20_micom_gene_to_genomic_region_v2.tsv', sep = '\t')
genomic.mapping <- genomic.mapping[!is.na(genomic.mapping$bp.start),]
genomic.mapping <- genomic.mapping[!duplicated(genomic.mapping$gene), ]



bool = genes %in% genomic.mapping$gene  #| genes %in% grch37.gnomad.2$gene #| genes %in% table.gnomad.2$gene[!is.na(table.gnomad.2$bp.end.gnomad.2)]
sum(bool)/length(bool) # ~ 0.58% not mapped\

#unmapped <- unique(genes[!bool])
#unmapped = na.omit(unmapped[unmapped %in% hgnc.aliases$to])


#write.table(unmapped, '~/Desktop/02JUN20_unmapped_genes.tsv', quote = F, row.names = F, col.names = F)



#colnames(table.ucsc)
#colnames(table.manuel)
#colnames(table.ncbi)


dfs = lapply(paths, function(path){
  
  write(paste0('mapping genomic regions ', path,'...'),stdout())
  
  df = read.csv(path, sep = '\t')
  df$gene = df$gene.first
  n1 = nrow(df)
  
  # merge columns and chrom
  df = merge(df, genomic.mapping, by.x = 'gene.final', by.y = 'gene', all.x = TRUE)
  n2 = nrow(df)
  
  if (n2 != n1) stop('expected all rows to be kept!')
  
  #newname
  #newpath = gsub('\\.tsv','\\.mapped\\.genomic\\.tsv', path)
  newpath = gsub('run031', 'run032', path)
  write.table(df, file = newpath, sep = '\t', quote = F, row.names = F)
  
  return(df)
})



if (F){

## get tables for referenc data.bases
grch37.ucsc <- read.csv('~/Desktop/GRCH37_gene_mapping_21MAY20.txt', sep = '\t')
grch37.ucsc.1 <- read.csv('~/Desktop/mart_export (2).txt', sep = '\t')
grch37.ucsc <- read.csv('data/01JUN20_ensemble_mart_export_genes.txt', sep = '\t')
grch37.refseq <- read.csv('~/Desktop/RefSeqMapped.tsv', sep = '\t')
grch37.gnomad <- read.csv('~/Desktop/gnomad/gnomad.v2.1.1.lof_metrics.by_transcript.txt', sep = '\t')
colnames(grch37.manuel.2) <- c('gene.from', 'gene.to', 'chrom', 'bp.start', 'bp.end')

# get table for ucsc
sum(genes %in% grch37.ucsc$UniProtKB.Gene.Name.ID)
table.ucsc = grch37.ucsc[match(genes, grch37.ucsc$UniProtKB.Gene.Name.ID), ]
table.ucsc = table.ucsc[complete.cases(table.ucsc), c(5,6,2,3)]
colnames(table.ucsc) <- c('gene', 'chrom.ucsc', 'bp.start.ucsc','bp.end.ucsc')

# get table for ucsc (refseq)
sum(genes %in% grch37.refseq$name2)
table.refseq = grch37.refseq[match(genes, grch37.refseq$name2), ]
table.refseq  = table.refseq[complete.cases(table.refseq), c(13,3,5,6)]#
table.refseq$chrom = gsub('chr','', table.refseq$chrom)
colnames(table.refseq)[1:4] <- c('gene', 'chrom', 'bp.start','bp.end')

# table for gnomad
sum(genes %in% grch37.gnomad$gene)
table.gnomad = grch37.gnomad[match(genes, grch37.gnomad$gene), ]
table.gnomad  = grch37.gnomad[,c(1, 76, 77, 78)]
colnames(table.gnomad)[1:4] <- c('gene', 'chrom.gnomad', 'bp.start.gnomad', 'bp.end.gnomad')
table.gnomad = table.gnomad[!duplicated(table.gnomad),]

# investigate difference in BP between refseq and UCSC
both = genes[genes %in% table.refseq$gene & genes %in% table.ucsc$gene]
head(table.ucsc[match(table.ucsc$gene, table.refseq$gene),])
mymerge = merge(table.ucsc, table.refseq, by = c('gene'))
diff1 = (mymerge$bp.start.refseq-mymerge$bp.start.ucsc)
diff2 = (mymerge$bp.end.refseq-mymerge$bp.end.ucsc)
boxplot(list(diff1, diff2), ylim = c(-2500,2500), ylab = 'BP difference')


### some genes are not found:

# Manuel data mining of gnomad (using geckodriver in python)
grch37.manuel <- readLines('~/Toolbox/dataminining/schema-webscraping/gnomad.tsv')
grch37.manuel <- lapply(grch37.manuel[2:length(grch37.manuel)], function(x) unlist(strsplit(x, '\t')))
hits = unlist(lapply(grch37.manuel, function(x){toupper(x[1]) == toupper(unlist(strsplit(x[2],'\\ ')))[1]}))
grch37.manuel <- grch37.manuel[hits]

# handle differential length 
table(unlist(lapply(grch37.manuel, function(x) length(x))))
sixes = unlist(lapply(grch37.manuel, function(x) length(x))) == 6
table.len.seven <- as.data.frame(do.call(rbind, grch37.manuel[!sixes]), stringsAsFactors = F)
table.len.six <- as.data.frame(do.call(rbind, grch37.manuel[sixes]), stringsAsFactors = F)
head(table.len.seven); head(table.len.six)

#
table.manuel <- data.frame(gene = c(table.len.seven$V1, table.len.six$V1), description =  c(table.len.seven$V2, table.len.six$V2),
                           gene.position = c(table.len.seven$V6, table.len.six$V5), stringsAsFactors = F)

# get chromosome position
table.manuel$gene.position <- gsub('Region ', '', table.manuel$gene.position)
table.manuel.regions <- as.data.frame(do.call(rbind, strsplit(table.manuel$gene.position, split = '\\:|\\-')), stringsAsFactors = F)
table.manuel <- cbind(table.manuel, table.manuel.regions)
colnames(table.manuel) <- c('gene', 'description', 'position','chrom.manuel', 'bp.start.manuel','bp.end.manuel')
write.table(table.manuel, '23MAY20_gnomad_table_automatically_curated.tsv', quote = F, row.names = F)
table.manuel$description <- NULL
table.manuel$position <- NULL


## manuel mapping of the remaining genes
grch37.ncbi <- read.csv('~/Desktop/manuel_mapping.csv', stringsAsFactors = F)
grch37.ncbi <- grch37.ncbi[!is.na(grch37.ncbi$Region),]
regions = lapply(strsplit(as.character(grch37.ncbi$Region), '\\ |\\, '), function(x) x[2])
regions = strsplit(gsub('\\(|\\)','', regions), split = '\\.\\.')
regions = as.data.frame(do.call(rbind, regions), stringsAsFactors = F)
colnames(regions) <- c('bp.start.ncbi','bp.end.ncbi')
grch37.ncbi = as.data.frame(cbind(grch37.ncbi, regions), stringsAsFactors = F)
grch37.ncbi$Region <- NULL
colnames(grch37.ncbi)[1:3] <- c('gene','chrom.ncbi', 'gene.new')
table.ncbi = grch37.ncbi

# deal with data-mined NCBI table 2
grch37.ncbi.2 <- read.csv('data/01JUN20_ncbi_webmining_gene_positions.tsv', sep = '\t')
grch37.ncbi.2$description <- NULL
grch37.ncbi.2$species <- NULL
grch37.ncbi.2$gene.found <- NULL

## deal with data-mined gnomad 2 table
grch37.gnomad.2 <- read.csv('data/01JUN20_gnomad_webmining_gene_positions.tsv', sep = '\t')
table.gnomad.2 <- grch37.gnomad.2
table.gnomad.2$chromosome.build <- NULL
table.gnomad.2$transcript.id <- NULL
table.gnomad.2$description <- NULL
table.gnomad.2$gene.id <- NULL
table.gnomad.2$reference <- NULL
table.gnomad.2$region <- NULL
table.gnomad.2$gene.found <- NULL
colnames(table.gnomad.2) <- c('gene','chrom.gnomad.2','bp.start.gnomad.2', 'bp.end.gnomad.2')

## deal with manuel mapping
grch37.manuel.2 <- read.table('data/01JUN20_final_mapped_genes.tsv', header = F)
table.manuel.2 <- grch37.manuel.2
colnames(table.manuel.2) <- c('gene', 'gene.to.manuel', 'chrom', 'bp.start', 'bp.end')



bool = genes %in% table.ucsc$gene | genes %in% table.gnomad.2$gene | genes %in% table.manuel$gene | 
  genes %in% table.ncbi$gene | genes %in% table.manuel.2$gene | genes %in% table.

bool = genes %in% genomic.mapping$gene

sum(bool)/length(bool)

genes[!bool]



####
# manuel 2


grch37.gnomad.2 <- read.csv('~/Toolbox/dataminining/schema-webscraping/02JUN20_gnomad.tsv', sep = '\t', header = F, na.strings=c(""," ","NA"))
grch37.gnomad.2.a <- grch37.gnomad.2[!grepl('References Ensembl', grch37.gnomad.2$V6), ]
colnames(grch37.gnomad.2.a) <- c('gene', 'description','chromosome.build', 'gene.id', 'transcript.id', 'region','reference')
grch37.gnomad.2.a <- grch37.gnomad.2.a[complete.cases(grch37.gnomad.2.a), ]
grch37.gnomad.2.b <- grch37.gnomad.2[grepl('References Ensembl', grch37.gnomad.2$V6), ]
colnames(grch37.gnomad.2.b) <- c('gene','chromosome.build', 'gene.id', 'transcript.id', 'region','reference','V7')
grch37.gnomad.2.b <- grch37.gnomad.2.b[complete.cases(grch37.gnomad.2.b), ]
grch37.gnomad.2.b$V7 <- NULL

grch37.gnomad.2 = merge(grch37.gnomad.2.a, grch37.gnomad.2.b, all = T)
grch37.gnomad.2$region <- gsub('Region ','', grch37.gnomad.2$region)

regions <- as.data.frame(do.call(rbind, strsplit(grch37.gnomad.2$region, split = '\\:|\\-')))
colnames(regions) <- c('chrom','bp.start','bp.end')

grch37.gnomad.2 <- cbind(grch37.gnomad.2, regions)

grch37.gnomad.2$gene.found <- unlist(lapply(strsplit(as.character(grch37.gnomad.2$description), '\\ '), function(x) x[1]))

grch37.gnomad.2[, sort(colnames(grch37.gnomad.2))]

bool = toupper(as.character(grch37.gnomad.2$gene)) == toupper(as.character(grch37.gnomad.2$gene.found))

# final
grch37.gnomad.2 <- grch37.gnomad.2[bool, ]
bool = genes %in% genomic.mapping$gene | genes %in% grch37.gnomad.2$gene | genes %in% dat$gene
sum(bool)/length(bool)

# get genes from NCBI

x0 <- genomic.mapping
x0 = as.data.frame(sapply(x0, as.character))


x1 <- grch37.gnomad.2[,c('gene','chrom','bp.start','bp.end')]
x1$genomic.location.origin <- 'GNOMAD'
x1$gene.approved <- NA
x1 = as.data.frame(sapply(x1, as.character))


x2 <- dat[,c('gene','chrom','bp.start','bp.end')]
x2$genomic.location.origin <- 'NCBI'
x2$gene.approved <- NA
x2 = as.data.frame(sapply(x2, as.character))

final <- merge(x0, x1, all = T)
final <- merge(final, x2, all = T)

final <- rbind(x0, x1, x2)

write.table(final, 'data/02JUN20_micom_gene_to_genomic_region_v2.tsv', sep = '\t',quote = F, row.names = F)


#x1 <- grch37.gnomad.2[complete.cases(grch37.gnomad.2),]

#x1$bp.start <- as.integer(x1$bp.start)
#x1$bp.end <- as.integer(x1$bp.end)
#x1$gene.found <- as.factor(x1$gene.found)
#x1$region <- as.factor(x1$region)

#x2 <- read.csv('data/01JUN20_gnomad_webmining_gene_positions.tsv', sep = '\t')


#x3[order(x3$gene),]
#x3 = rbind(x1, x2)
#x3 = x3[!duplicated(x3),]


#write.table(x3, '01JUN20_gnomad_webmining_gene_positions.tsv', row.names = F, quote = F, sep = '\t')

## NCBI

#a <- read.csv('~/Toolbox/dataminining/schema-webscraping/02JUN20_ncbi_mapped_genes_a.tsv', sep = '\t', header = T, na.strings=c(""," ","NA"))
#b <- read.csv('~/Toolbox/dataminining/schema-webscraping/02JUN20_ncbi_mapped_genes_b.tsv', sep = '\t', header = T, na.strings=c(""," ","NA"))

#dat = rbind(a, b)
#dat <- dat[complete.cases(dat),]
#dat$gene.found <- unlist(lapply(strsplit(as.character(dat$gene), '\\ '), function(x) x[1]))
#colnames(dat)[1:2] <- c('gene', 'description')

#dat = dat[complete.cases(dat),]
#dat$gene == dat$gene.found

#dat[,c('gene','chrom','bp.start','bp.end')]




#write.table(dat, '01JUN20_ncbi_webmining_gene_positions.tsv', row.names = F, quote = F, sep = '\t')


#a <- read.csv('~/Toolbox/dataminining/schema-webscraping/01JUN20_ncbi_mapped_genes_a.tsv', sep = '\t', header = T, na.strings=c(""," ","NA"))
#b <- read.csv('~/Toolbox/dataminining/schema-webscraping/01JUN20_ncbi_mapped_genes_b.tsv', sep = '\t', header = T, na.strings=c(""," ","NA"))
#c <- read.csv('~/Toolbox/dataminining/schema-webscraping/01JUN20_ncbi_mapped_genes_c.tsv', sep = '\t', header = T, na.strings=c(""," ","NA"))

#d <- rbind(a, b, c)
#d <- d[complete.cases(d),]
#d$gene.found <- unlist(lapply(strsplit(as.character(d$gene), '\\ '), function(x) x[1]))
#d[d$gene_query != d$gene.found,]
#colnames(d) <- c('gene.query','description','chrom','bp.start','bp.end','species','gene.found')
#write.table(d, '01JUN20_ncbi_webmining_gene_positions.tsv', row.names = F, quote = F, sep = '\t')

}





