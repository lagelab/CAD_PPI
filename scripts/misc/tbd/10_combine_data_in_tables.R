options(java.parameters = "-Xmx8000m")

library(xlsx)

# get the files that passed QC metrics
accepted = readLines('~/Projects/03_MICOM/21JUL20_micom_qced_ip_paths.txt')
files = list.files('~/Projects/03_MICOM/data/genoppi_input/run019/')
files.path.full = list.files('~/Projects/03_MICOM/data/genoppi_input/run019/', full.names = T)
target = files.path.full[files %in% accepted]

if (length(target) != 18) stop('expected 18 datasets!')

# load all IPs 
ips = lapply(target, read.csv, sep = '\t')
names(ips) = basename(target)

# prep for naming
path_bait_table <- read.csv('data/pathname_bait_to_gene.csv', stringsAsFactors = F)
translate_path_bait <- function(x, table) {
  return(table$To[table$From == x])
}

# get new names
newnames <- unlist(lapply(names(ips), function(path){
  path.bait.direct = unlist(lapply(strsplit(path,'\\.|(vs)'), function(x) x[2]))
  path.bait = translate_path_bait(path.bait.direct, path_bait_table)
  path.cell = unlist(lapply(strsplit(path, '\\.'), function(x) x[3]))
  outname = paste(path.cell,'.',path.bait,'vsControl', sep = '')
  return(outname)
}))
  
# rename
oldnames <- names(ips)
names(ips) <- newnames
names(ips)[9] <- paste0(names(ips)[9], '(I)') 
names(ips)[10] <- paste0(names(ips)[10], '(I)')
ordered.names <- names(ips)[order(names(ips))]

outfile = 'Supplementary_table_1b.xlsx'
for (outname in ordered.names){
  print(outname)
  write.xlsx(ips[[outname]], file=outfile, sheetName=outname, append=TRUE, row.names=FALSE)
}

newnames <- unlist(lapply(names(ips), function(path){
  path.bait.direct = unlist(lapply(strsplit(path,'\\.|(vs)'), function(x) x[2]))
  path.bait = translate_path_bait(path.bait.direct, path_bait_table)
  path.cell = unlist(lapply(strsplit(path, '\\.'), function(x) x[3]))
  outname = paste(path.cell,'.',path.bait,'vsControl', sep = '')
}))

# make description
ordered.oldnames <- oldnames[order(names(ips))]

description.names <- unlist(lapply(ordered.oldnames, function(path){
  
  path.bait.direct = unlist(lapply(strsplit(path,'\\.|(vs)'), function(x) x[2]))
  path.bait = translate_path_bait(path.bait.direct, path_bait_table)
  path.bait.inweb <- gsub('PLPP3', 'PPAP2B', gsub('JCAD', 'KIAA1462',path.bait))
  path.vs = unlist(lapply(strsplit(path,'\\.|(vs)'), function(x) x[3]))
  path.vs[path.vs == 'Whitehead'] <- NA
  path.cell = unlist(lapply(strsplit(path, '\\.'), function(x) x[3]))
  path.date = unlist(lapply(strsplit(path, '\\.'), function(x) x[5]))
  path.facility = ifelse(grepl('Broad', path), 'Broad', 'Whitehead')
  path.type = unlist(lapply(strsplit(path, '\\.'), function(x) x[4]))
  outname = paste0(path.bait,' vs Control in ', path.cell, ' (', path.type,')')
  return(outname)
  
}))







