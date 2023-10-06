setwd('~/Projects/03_MICOM/')
# load paths
#devtools::load_all('~/Projects/08_genoppi_lagelab/Genoppi/')

# Get breakdown of what was removed and why
paths <- list.files('data/genoppi_input/run018', full.names = T)
run = 'run018'

breakdown = lapply(paths, function(path) {
  
  # get raw metrics from data
  df = read.csv(path, sep = '\t', stringsAsFactors = F)
  
  # get data/path origin
  path.bait.direct = unlist(lapply(strsplit(path,'\\.|(vs)'), function(x) x[2]))
  path.bait = translate_path_bait(path.bait.direct, path_bait_table)
  path.vs = unlist(lapply(strsplit(path,'\\.|(vs)'), function(x) x[3]))
  path.vs[path.vs == 'Whitehead'] <- NA
  path.cell = unlist(lapply(strsplit(path, '\\.'), function(x) x[3]))
  path.date = unlist(lapply(strsplit(path, '\\.'), function(x) x[5]))
  path.facility = ifelse(grepl('Broad', path), 'Broad', 'Whitehead')
  
  # tally data
  rows.total <- nrow(df)
  rows.discarded <- sum(df$remove)
  what = unlist(strsplit(df$rmflag, split = '\\;'))
  what = what[what != '']
  dat <- cbind(run = run, path = basename(path), facility = path.facility, bait = path.bait, vs = path.vs, cell = path.cell, data.frame(table(what), stringsAsFactors = F))
  
  # get summary stats
  dat$rows.total = rows.total
  dat$rows.discarded = rows.discarded
  dat$pct.of.total = round(dat$Freq / dat$rows.total * 100, 1)
  dat$pct.of.discarded = round(dat$Freq / dat$rows.discarded * 100, 1)
  colnames(dat) <- c('run','path', 'facility','bait','vs','cell','reason', 'freq', 'rows (total)', 'rows discarded', ' % of total', '% of discarded')
  dat <- rbind(dat, '')
  
  return(dat)
})


newdat <- do.call(rbind, breakdown)
write.csv(newdat, '~/Desktop/21JUL20_MICOM_run019_QC_SUMMARY.csv', row.names = F)


# Find distribution of ribosomal proteins
paths <- list.files('data/genoppi_input/run019', full.names = T)
ribo <- lapply(paths, function(path) {
  df = read.csv(path, sep = '\t', stringsAsFactors = F)
  df = df[df$significant, ]
  return(sum(grepl('^RPL', df$gene)/nrow(df)))
})
names(ribo) <- paths
hist(unlist(ribo), xlim = c(0, 1), breaks = 15, main = '% Ribosomal proteins in interactors', xlab = '% of enriched proteins')
abline(v = median(unlist(ribo)), col = 'red')





# check some data

path = 'data/genoppi_input/run017/Genoppi_Whitehead_martin6_MinImputed.HDAC9vsMock.SMC.3xFlag.20JUL2020.tsv'

df = read.csv(path, sep = '\t') %>% id_enriched_proteins()

df %>% plot_scatter_basic() 

# look up weird entry

path = 'data/genoppi_input/run017/Genoppi_Whitehead_martinBroad_MinImputed.PHACTR1vsMock.SMC.expanded.20JUL2020.tsv'

df = read.csv(path, sep = '\t') %>% id_enriched_proteins()


