setwd('~/Projects/03_MICOM/')

#
paths = list.files('data/genoppi_input/run032', full.names = T)
paths.whitehead = !grepl('Broad', paths)
paths.broad = grepl('Broad', paths)
hgnc.remap <- read.csv('data/02JUN20_hgnc_hgnc_prevsymbol_to_current_symbol.tsv', sep = '\t', stringsAsFactors = F)
baits = c('BCAS3' ,'FLT1', 'KCNK5', 'ARHGEF26', 'JCAD', 'FN1', 'EDNRA', 'HDAC9', 'PLPP3', 'ADAMTS7','PHACTR1', 'KSR2')

# read lists of common contanimants
#contaminants <- readLines('data/gi_contaminants.fasta')
#contaminants <- contaminants[grepl('(gi)|(sp)',contaminants)]

#known.contaminants = read.csv('data/17JUN20_thegpm_crap.tsv', sep ='\t')
#known.contaminants = known.contaminants[known.contaminants$table == 2, ]
#known.contaminants$species = unlist(lapply(strsplit(as.character(known.contaminants$id), '\\_'), function(x) x[2]))
#known.contaminants = known.contaminants[known.contaminants$species == 'HUMAN',]

#dprot <- read.delim(url('https://www.uniprot.org/mapping/M20200617A94466D2655679D1FD8953E075198DA800804FS.tab'))
#dprot = merge(df, dprot, by.x = 'id', by.y = 'From')
#colnames(dprot) <- c('id','table', 'description', 'reason', 'gene')
#write.table(dprot, 'data/17JUN20_thegpm_crap.tsv', sep ='\t', quote = F, row.names = F)




is_human_or_na <- function(x){
  bool = rep(F, length(x))
  b1 = bool
  b2 = bool
  b1[is.na(x)] <- T
  b2[x == 'HUMAN'] <- T 
  return(b1 | b2)
}

# refine columns and remove redundant ones
dfs = lapply(paths, function(path){
  
  #print(path)
  df = read.csv(path, sep = '\t', stringsAsFactors = F)
  #print(head(df))
  
  path.bait = unlist(lapply(strsplit(path,'\\.|(vs)'), function(x) x[2]))
  path.bait = gsub('(\\(WT\\))|(\\(MT\\))|(\\-beta)|(\\-alpha)', '', path.bait)
  
  #--------------------------------------
  # standardize column names
  
  if ('Accession' %in% colnames(df)){
    df$accession <- df$Accession
  } else if ('accession_number' %in% colnames(df)){
    df$accession <- df$accession_number
  } else {
    df$accession <- NA
  }
  
  if ('ec.rep1' %in% colnames(df)){
    df$rep1 <- df$ec.rep1
    df$rep2 <- df$ec.rep2
  }
  
  if ('smc.rep1' %in% colnames(df)){
    df$rep1 <- df$smc.rep1
    df$rep2 <- df$smc.rep2
  }
  
  if (!'imputed' %in% colnames(df)){
    df$imputed <- NA
  }
  
  if (!'flag.tag' %in% colnames(df)){
    df$flag.tag <- NA
  }
  
  if (!'flag.tag' %in% colnames(df)){
    df$flag.tag <- NA
  }
  
  #---------------------------------------
  # add tag to remove proteins that are 
  # * uncharcaterized
  # * not human
  # * predicted to not be human
  # * trypsins used in MS digestion
  # * rep columns that are NAs
  # * not mapped to a chromosomal location
  # * gene names are NA
  # * keratins (Regex by KRT)
  
  
  df$ignore <- ifelse('flag.tag' %in% colnames(df), !is.na(df$flag.tag), F) | df$gene.final %in% path.bait
  df$remove <- FALSE
  df$rmflag <- ''
  
  
  if ('Uncharacterized' %in% colnames(df)){
    
    #if (any(df$Uncharacterize)) browser()
    #df = df[!df$Uncharacterized | df$ignore, ]
    df$remove <- df$remove | df$Uncharacterized
    df$rmflag[df$Uncharacterized] <- paste(df$rmflag[df$Uncharacterized], 'Uncharacterized', sep = ';')
  }
  
  if ('species' %in% colnames(df)){
    
    #df = df[df$species == 'HUMAN' | is.na(df$species) | df$ignore,]
    df$remove <- df$remove | !is_human_or_na(df$species)
    df$rmflag[!is_human_or_na(df$species)] <- paste(df$rmflag[!is_human_or_na(df$species)], 'NA species', sep = ';')
    
  }
  
  if ('Description' %in% colnames(df)){
    bool.animal <- grepl('(Porcine)|(Capra hircus)|(Bos taurus)', df$Description)
    #df = df[!bool.animal | df$ignore, ]
    df$remove <- df$remove | bool.animal
    df$rmflag[bool.animal] <- paste(df$rmflag[bool.animal], 'non-HUMAN', sep = ';')
  }
  
  if ('X.Unique' %in% colnames(df)){
    
    #df = df[df$X.Unique > 1 | df$ignore, ]
    df$unique.peptides <- df$X.Unique
    df$remove <- df$remove | df$unique.peptides <= 1
    
  } else if ('unique.peptides' %in% colnames(df)){
    
    #df = df[df$unique.peptides > 1 | df$ignore, ]
    df$remove <- df$remove | df$unique.peptides <= 1
    
  } else if ('ec.unique.peptides' %in% colnames(df)){
    
    #df = df[df$ec.unique.peptides > 1 | df$ignore, ]
    df$unique.peptides <- df$ec.unique.peptides
    df$remove <- df$remove | df$unique.peptides <= 1
    
  } else if ('smc.unique.peptides' %in% colnames(df)){
    
    #df = df[df$smc.unique.peptides > 1 | df$ignore, ]
    df$unique.peptides <- df$smc.unique.peptides
    df$remove <- df$remove | df$unique.peptides <= 1
    
  } else {
    df$unique.peptides <- NA
  }
  
  if ('unique.peptides' %in% colnames(df)){
    df$rmflag[df$unique.peptides <= 1] <- paste(df$rmflag[df$unique.peptides <= 1 & !is.na(df$unique.peptides)], 'unique peptides < 2', sep = ';')
  }
  
  rep.na <- is.na(df$rep1) | is.na(df$rep2)
  if (any(rep.na)){
    df$remove <- df$remove | rep.na
    df$rmflag[rep.na] <- paste(df$rmflag[rep.na], 'replicates are NA', sep = ';')
  }
  
  chrom.na <- is.na(df$chrom)
  if (any(chrom.na)){
    df$remove <- df$remove | chrom.na
    df$rmflag[chrom.na] <- paste(df$rmflag[chrom.na], 'chromosomes are NA', sep = ';')
  }
  
  gene.na <- is.na(df$gene.final)
  if (any(gene.na)){
    df$remove <- df$remove | gene.na
    df$rmflag[gene.na] <- paste(df$rmflag[gene.na], 'genes are NA', sep = ';')
  }
  
  keratins <- grepl('^KRT[0-9]', df$gene.final) | grepl('^KRTAP', df$gene.final)
  if (any(keratins)){
    df$remove <- df$remove | keratins
    df$rmflag[keratins] <- paste(df$rmflag[keratins], 'Keratin or keratin-associated protein', sep = ';')
  }
  
  # remove MS duplicated entries, e.g. unresolved isoforms.
  dups <- duplicated(df[,c('gene.final','rep1', 'rep2')])
  if (any(dups)){
    df$remove <- df$remove | dups
    df$rmflag[dups] <- paste(df$rmflag[dups], 'unresolved isoform', sep = ';')
  }
  
  
  #contaminants <- df$gene.final %in% known.contaminants$gene
  #if (any(contaminants)){
  #  df$remove <- df$remove | contaminants
  #}
  
  
  df$remove <- df$remove & !df$ignore
  df$rmflag[!df$remove] <- ''
  
  
  #--------------------------------------
  # some genenames are old (not approved)
  # let's remap them.
  
  na.bool = is.na(df$gene.approved)
  if (any(!na.bool)){
    
    # deal with unusual cases
    df$gene.approved[!na.bool] = gsub('SPINDOC\xca', 'SPINDOC', df$gene.approved[!na.bool])
    
    old = df$gene.final[!na.bool]
    new = df$gene.approved[!na.bool]
    
    # some qc metrics
    cat('[START] re-naming gene.final ..\n')
    print(data.frame(from=old, to=new))
    cat('[END]\n')
    
    df$gene.final[!na.bool] <-df$gene.approved[!na.bool]
    
  }
  
  hgnc.bool = df$gene.final %in% hgnc.remap$from
  # hgnc remap to approved symbol
  if (any(hgnc.bool)){
    
    from <- as.character(df$gene.final[hgnc.bool])
    to <- unlist(lapply(from, function(x) as.character(hgnc.remap$to[hgnc.remap$from == x])))
    df$gene.final[hgnc.bool] <- to
    
  }
  
  # newname
  newpath = gsub('\\.tsv', paste0('\\.',toupper(format(Sys.time(), "%d%b%Y")),'\\.tsv'), path)
  newpath = gsub('run032', 'run033', newpath)
  write.table(df, file = newpath, sep = '\t', quote = F, row.names = F)
  
  print(head(df))
  
  return(df)
})






