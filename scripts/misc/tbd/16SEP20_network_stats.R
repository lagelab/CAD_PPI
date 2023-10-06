

combine_int <- function(f){
  df = read.csv(f, sep = '\t')
  df = df[ ,c('Bait', 'GeneName', 'Interactor', 'CellType')]
  colnames(df) <- c('Bait','Gene','Significant','Cell')
  return(df)
}

# table for EC
path = 'derived/maps/maps-run019/'
filesEC <- list.files(path,'EC\\.InteractorTable\\.txt', full.names = T)
filesEC <- filesEC[!grepl('COMBINED',filesEC)]
tabEC <- do.call(rbind, lapply(filesEC, combine_int))

## Table for SMC
path = 'derived/maps/maps-run019/'
filesSMC <- list.files(path,'\\_SMC\\.InteractorTable\\.txt', full.names = T)
filesSMC <- filesSMC[!grepl('COMBINED',filesSMC)]
tabSMC<- do.call(rbind, lapply(filesSMC, combine_int))
tabSMC$Bait[tabSMC$Bait == 'JCAD'] <- 'KIAA1462' # JCAD is not found in inweb..

## Table for SMC
path = 'derived/maps/maps-run019/'
filesECSMC <- list.files(path,'EC\\-SMC\\.InteractorTable\\.txt', full.names = T)
filesECSMC <- filesECSMC[!grepl('COMBINED',filesECSMC)]
tabECSMC<- do.call(rbind, lapply(filesECSMC, combine_int))

## how many interactions are found in InWeb
tabs <- list(tabEC, tabSMC, tabECSMC)
names(tabs) <- c("EC", "SMC", 'EC-SMC')

# check baits
checkinweb = unlist(lapply(unique(c(tabSMC$Bait, tabEC$Bait)), function(bait) any(get_inweb_list(bait)$significant)))
if (any(!checkinweb)) stop('at least one bait did not exist in inweb!')

haarst <- read.table('genelist/27AUG20_haarst2017_protein_coding_genelist.tsv', sep = '\t', stringsAsFactors = F, header = T)


# what interactions are actually verfied?
interactions <- lapply(tabs, function(df){
  
  baits = unique(df$Bait)
  df$ininweb <- F
  for (b in baits){
    inweb = get_inweb_list(b)
    which = df$Bait %in% b & df$Gene %in% inweb$gene[inweb$significant]
    if (any(which)) df[which, ]$ininweb <- T
  }
  df$haarst <- df$Gene %in% haarst$geneName
  
  # only keep interactors
  df = df[df$Significant, ]
  # remove duplicates
  #df = df[!duplicated(df),]
  
  return(df)
  
})

# how many total interactions
lapply(interactions, function(x) nrow(x))

# how many inweb interactions
lapply(interactions, function(x) sum(x$ininweb))

# how many are haarst risk genes?
lapply(interactions, function(x) sum(x$haarst))

# get counts for graph interaction partners
table(as.vector(table(interactions$EC$Gene)))
table(as.vector(table(interactions$SMC$Gene)))
table(as.vector(table(interactions$`EC-SMC`$Gene)))

## table with combined interactions

combinedEC <- list.files(path,'EC\\.InteractorTable\\.txt', full.names = T)
combinedEC <- filesEC[grepl('COMBINED',combinedEC)]
cEC <-  read.csv(combinedEC, sep = '\t')

# interactions subsetted



