# Get Haarst 2017 table data
table.xx <- read.csv('data/credible genesets/312086_table_xx.csv', stringsAsFactors = F)
martin <- read.csv('data/credible genesets/15MAY20_Martins_causal_genes.csv', stringsAsFactors = F)
colnames(martin)[1] <- 'rsID'
martin$rsID <- gsub(' ', '', martin$rsID)
variants = unlist(strsplit(martin$rsID, split = '\\, '))

# get posterior probabilities
martin.c <- merge(martin, table.xx, by = 'rsID', all.x = T)
martin.c.subset <- martin.c[, c('LocusID','rsID', 'Posterior.probability.for.causaility', 'martin.casual.gene')]

## biomart queries
library(biomaRt)
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
attributes = listAttributes(ensembl)
mart = getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'strand', 'chromosome_name', 'start_position', 'end_position','ensembl_gene_id_version'), 
      filters = 'chromosome_name', 
      values = c(1:22,'X','Y'), 
      mart = ensembl)
colnames(mart) <- c('hgnc.symbol', 'Gene.stable.ID','strand.dir','chr.name','Gene.start..bp.', 'Gene.end..bp.', 'GENCODE_id')

# merge ensemble_gene_id with 
grch37.hgnc <- mart #<- merge(mart, grch37, by = 'Gene.stable.ID',)
table.xx2 = merge(table.xx, mart, all.x = T, by = 'GENCODE_id')

# variant effect predictor
vep <- read.table('data/credible genesets/28MAY20_MICOM_SNPs_VEP_2.txt', sep = '\t', header = T)
colnames(vep)[1] <- 'rsID'




# get locus and nearby 
get_locus <- function(id, width = 100000, table = table.xx2, gene.name.size = 2.5, martin.snp = '', martin.bait = '', title = ''){
  
  # require data 
  require(ggplot2)
  require(ggrepel)
  require(gtable)
  require(grid)
  require(RColorBrewer)
  
  # query from locus
  locus = table[table$LocusID %in% id, c(colnames(table)[1:15])]
  locus = locus[!is.na(locus$Chromosome), ]
  locus$martinsnp = as.factor(ifelse(locus$rsID %in% martin.snp, 'Y', 'N'))
  locus$label = locus$rsID
  
  if (nrow(locus) < 1) stop(paste('locus',id,'is not valid!'))
  
  # data VEP: combine CONSEQUENCE with SYMBOL
  data.vep = merge(locus, vep, all.x = T, by = 'rsID')
  consequence = unlist(lapply(unique(data.vep$rsID), function(id){ 
    
    cons.symbol = data.vep[data.vep$rsID %in% id, c('Consequence', 'SYMBOL')]
    cons.symbol$SYMBOL[cons.symbol$SYMBOL == '-'] <- NA
    cons.symbol = as.data.frame(cons.symbol[!duplicated(cons.symbol),], stringsAsFactors = F)
    cons.symbol = as.character(apply(cons.symbol, 1, paste, collapse = '='))
    
    return(paste(cons.symbol, collapse = '\n'))
    
    }))
  
  locus$consequence = consequence# combine with gene
  locus$label = paste(locus$rsID, '\nVEP:[', unique(locus$consequence) ,']')
  
  if (nrow(data.vep) < 1) stop(paste(' VEP data is not valid!'))
  
  # get limits
  genomic.position <- as.vector(na.omit(locus$MarkerName..chr.pos.alleleA_alleleB.))
  position <- as.numeric(unlist(lapply(strsplit(genomic.position, '\\:|\\_'), function(x) x[2])))
  chrom <- unique(as.numeric(unlist(lapply(strsplit(genomic.position, '\\:|\\_'), function(x) x[1]))))
  lims = list(min = min(position) - width, max = max(position) + width)
  
  #Create a custom color scale
  cols <- brewer.pal(5,"Set1")[1:2]
  names(cols) <- c('Y','N')
  colScale <- scale_colour_manual(name = "is.bait",values = cols)
  
  
  # setup plotting for SNPs
  p1 <- ggplot(locus, aes(x=Position..SNP.in.LD, y=Posterior.probability.for.causaility, color = martinsnp)) +
    geom_point() + 
    scale_x_continuous(limits = c(lims$min,lims$max)) +
    xlab('') + 
    ylab('Posterior probability of causaility') +
    ggtitle(paste(paste('locus',id), title), paste0('Lead SNP (',locus$Gwas.SNP.1[locus$Gwas.SNP],') LD range = [', round(min(locus$R2),2), ', ', round(max(locus$R2), 2),']')) +
    theme_bw() + 
    guides(color=guide_legend(title="Is MICOM SNP?"))
  
  # map variant names and consequence
  p1 <- p1 + geom_text_repel(aes(x=Position..SNP.in.LD, y=Posterior.probability.for.causaility, label=label),size = 2)
  p1 <- p1 + colScale
  
  
  # get genes within limits
  genes = grch37.hgnc[((grch37.hgnc$Gene.end..bp. >= lims$min & grch37.hgnc$Gene.end..bp. < lims$min) |
              (grch37.hgnc$Gene.end..bp. >= lims$min & grch37.hgnc$Gene.end..bp. < lims$max) |
              (grch37.hgnc$Gene.start..bp. <= lims$max  & grch37.hgnc$Gene.start..bp. > lims$min)) & 
                grch37.hgnc$chr.name %in% as.character(chrom), 
              c('hgnc.symbol', 'chr.name','Gene.start..bp.', 'Gene.end..bp.','strand.dir')]
  genes = genes[!duplicated(genes) & genes$hgnc.symbol != '',]
  
  if (nrow(genes) > 0){
  
    # order genes
    genes = genes[order(genes$Gene.start..bp.), ]
    genes$ycenter = 1:nrow(genes)
    genes$ymin = genes$ycenter-0.5
    genes$ymax = genes$ycenter+0.5
    genes$len = genes$Gene.end..bp.-genes$Gene.start..bp.
    genes$label = paste(genes$hgnc.symbol, ifelse(genes$strand.dir == 1, '(+)', '(-)' ))
    genes$label = paste(genes$label, '[', genes$len,'bases]')
    genes$is.bait = as.factor(ifelse(genes$hgnc.symbol %in% martin.bait, 'Y', 'N'))
    genes$gene.start.bp.cut = unlist(lapply(genes$Gene.start..bp., function(x) max(x, lims$min)))
    genes$gene.end.bp.cut = unlist(lapply(genes$Gene.end..bp., function(x) min(x, lims$max)))
    genes$xcenter = (genes$gene.start.bp.cut+genes$gene.end.bp.cut )/2

    # plot gene position
    p2 <- ggplot(genes, aes(xmin=gene.start.bp.cut, xmax = gene.end.bp.cut, ymin = ymin, ymax = ymax, x = genes$xcenter, y = genes$ycenter, label = label, fill = is.bait)) +
      geom_rect(alpha=0.2, color = 'black') +
      geom_text(size = gene.name.size) +
      scale_x_continuous(limits = c(lims$min,lims$max)) +
      xlab(paste0('BP Position (chromosome ',chrom,')')) +
      ylab(paste('Genes, n =',nrow(genes))) + 
      guides(fill=guide_legend(title="Is MICOM bait?")) +
      theme_bw() +
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    
    p2 <- p2 + colScale
  
  } else {
    warning(paste('widths are not big enough for ID', id))
    p2 = NULL
  }
  

  return(list(probs=p1, genes=p2, genes = genes))
}

draw_stacked_plot <- function(loc){
  
  require(grid)
  require(gridExtra)
  require(cowplot)
  
  grid.newpage()
  if (!is.null(loc[[1]]) & !is.null(loc[[2]])){
    algn = align_plots(loc[[1]], loc[[2]], align="hv", axis="tblr")
    grid.arrange(algn[[1]], algn[[2]], ncol = 1)
    return(invisible(T))
  } else {
    print(loc[[1]])
    return(invisible(F))
  }
}


bait = c(martin$martin.casual.gene, 'KIAA1462', 'PPAP2B')

pdf('28MAY20_micom_loci_summaries_01.pdf', width = 9, height = 12)

loc = get_locus(1286, width = 50000, martin.snp = martin$rsID, martin.bait = bait, title = 'Suggested causal: KCNK5')
draw_stacked_plot(loc)

# draw plot (KCNK5)
loc = get_locus(976, width = 50000, martin.snp = martin$rsID, martin.bait = bait, title = 'Suggested causal: KCNK5')
draw_stacked_plot(loc)

# draw plot (FN1)
loc = get_locus(1242, width = 50000, martin.snp = martin$rsID, martin.bait = bait, title = 'Suggested causal: FN1')
draw_stacked_plot(loc)

# draw plot (EDNRA)
loc = get_locus(1286, width = 50000, martin.snp = martin$rsID, martin.bait = bait, title = 'Suggested causal: EDNRA')
draw_stacked_plot(loc)

# draw plot (JCAD)
loc = get_locus(1211, width = 50000, martin.snp = martin$rsID, martin.bait = bait, title = 'Suggested causal: JCAD (KIAA1462)')
draw_stacked_plot(loc)

# draw plot (PPAP2B ~ PLPP3)
loc = get_locus(1354, width = 50000, martin.snp = martin$rsID, martin.bait = bait, title = 'Suggested causal: PLPP3 (PPAP2B)')
draw_stacked_plot(loc)

# draw plot (EDN1, PHACTR1)
loc = get_locus(1450, width = 1000000, martin.snp = martin$rsID, martin.bait = bait, title = 'Suggested causal: EDN1 and PHCATR1')
draw_stacked_plot(loc)

# draw plot (FLT)
loc = get_locus(934, width = 50000, martin.snp = martin$rsID, martin.bait = bait, title = 'FLT')
draw_stacked_plot(loc)

## SNP not found

# draw plot (ARHGEF26)
loc = get_locus(1091, width = 50000, martin.snp = martin$rsID, martin.bait = bait, title = 'ARHGEF26')
draw_stacked_plot(loc)

# draw plot (ADAMTS7)
loc = get_locus(1375, width = 50000, martin.snp = martin$rsID, martin.bait = bait, title = 'ADAMTS7')
draw_stacked_plot(loc)

# (HDAC9)
loc = get_locus(1305, width = 50000, martin.snp = martin$rsID, martin.bait = bait, title = 'HCAD9')
draw_stacked_plot(loc)

# (EDN) Note, this is only with 
loc = get_locus(562, width = 50000, martin.snp = martin$rsID, martin.bait = bait, title = 'HCAD9')
draw_stacked_plot(loc)

loc = get_locus(1450, width = 750000, martin.snp = martin$rsID, martin.bait = bait, title = 'HCAD9')
draw_stacked_plot(loc)




graphics.off()

lst = list()


## generate summary table
pdf('28MAY20_micom_loci_selected_summaries_500000bp.pdf', width = 9, height = 12)
for (id in na.omit(unique(table$LocusID))){
  #print(id)
  loc = get_locus(id, width = 500000, martin.snp = martin$rsID, martin.bait = bait, title = 'ID')
  if (any(loc$genes$data$is.bait == 'Y') | any(loc$probs$data$rsID %in% martin.c.subset$rsID)){
    draw_stacked_plot(loc)
    
    # get micom gene information 
    micom.gene = loc$genes$data$hgnc.symbol[loc$genes$data$is.bait == 'Y']
    vec <- c(vec, micom.gene)
    print(paste(id, '-',micom.gene))
    SNPs = paste(loc$probs$data$rsID, collapse = ', ')
    effect = paste(loc$probs$data$consequence, collapse = ', ')
    
    # get posterior probablity
    prob.posterior <- loc$probs$data$Posterior.probability.for.causaility
    most.probable.snp <- loc$probs$data$rsID[max(prob.posterior) == prob.posterior]
    most.probable.snp.posterior <- loc$probs$data$Posterior.probability.for.causaility[max(prob.posterior) == prob.posterior]
    most.probable.snp.consequence <- loc$probs$data$consequence[max(prob.posterior) == prob.posterior]
    most.probable.snp.location <- loc$probs$data$Position..SNP.in.LD[max(prob.posterior) == prob.posterior]
    
    # gene and closets gene
    nearby.genes = paste(loc$genes$data$hgnc.symbol, collapse = ', ')
    n.genes.nearby = length(loc$genes$data$hgnc.symbol)
    bp.diff = abs(most.probable.snp.location - loc$genes$data$xcenter)
    closets.gene = loc$genes$data$hgnc.symbol[min(bp.diff) == bp.diff]

    row = data.frame(id, SNPs, most.probable.snp, most.probable.snp.posterior, gsub('\n',', ',most.probable.snp.consequence), n.genes.nearby ,nearby.genes, closets.gene, micom.gene)
    lst[[id]] = row
  
  }
}
graphics.off()



lst[sapply(lst, is.null)] <- NULL
dat = do.call(rbind, lst)
dat = as.data.frame(dat)
colnames(dat) <- c('locus', 'SNPs', 'most.probable.snp',' most.probable.snp.posterior',' most.probable.snp.consequence', 'n.genes.nearby' ,'nearby.genes', 'closets.gene', 'micom.gene')

write.table(dat, file = 'micom_summary_snps_500000bp.tsv', quote = F, row.names = F, sep = '\t')














  
  