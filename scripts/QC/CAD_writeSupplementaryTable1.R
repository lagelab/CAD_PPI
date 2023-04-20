setwd('~/Projects/03_MICOM/')

# ready paths and prefixes
date = '201120'
run = 'run056'

# setup paths and extract tier 1 paths
prefix = paste0(date,'_',run,'_micom')
paths = list.files(paste0('data/genoppi_input/', run), full.names = T)
summary <- read.csv('201120_run056_micom_tier_summary.tsv', sep = '\t')
summary <- summary[!is.na(summary$tier),]

# filter out things we don't want to analyze
summary$ok <- 
  summary$tier == 1 &
  !summary$Bait %in% 'KCNK5' &
  !summary$Cell.line %in% 'HEK293T' &
  !grepl('\\(MT', summary$Data.path) &
  !grepl('\\.SGS', summary$Data.path)

# Datasets included in run19 (18 rows)
# Genoppi_Whitehead_martin10_MinImputed.FN1vsMock.SMC.Endogenous.21JUL2020.tsv				
# Genoppi_Whitehead_martin10_MinImputed.JCADvsMock.SMC.3xFLAG.21JUL2020.tsv				
# Genoppi_Whitehead_martin11_MinImputed.FN1vsMock.EC.Endogenous.21JUL2020.tsv				
# Genoppi_Whitehead_martin12_MinImputed.EDN1vsMock.EC.3xFLAG.21JUL2020.tsv				
# Genoppi_Whitehead_martin12_MinImputed.PHACTR1vsMock.SMC.3xFLAG.21JUL2020.tsv				
# Genoppi_Whitehead_martin13_MinImputed.PHACTR1vsMock.EC.3xFLAG.21JUL2020.tsv				
# Genoppi_Whitehead_martin13_MinImputed.PLPP3vsMock.SMC.3xFLAG.21JUL2020.tsv				
# Genoppi_Whitehead_martin14_MinImputed.EDNRAvsMock.SMC.3xFLAG.21JUL2020.tsv				
# Genoppi_Whitehead_martin15_MinImputed.ARHGEF26(WT)vsMock.SMC.3xFLAG.21JUL2020.tsv				
# Genoppi_Whitehead_martin15_MinImputed.FLT1vsMock.EC.3xFLAG.21JUL2020.tsv				
# Genoppi_Whitehead_martin15_MinImputed.PLPP3vsMock.EC.3xFLAG.21JUL2020.tsv				
# Genoppi_Whitehead_martin5_MinImputed.BCAS3vsMock.SMC.3xFlag.21JUL2020.tsv				
# Genoppi_Whitehead_martin6_MinImputed.HDAC9vsMock.SMC.3xFlag.21JUL2020.tsv				
# Genoppi_Whitehead_martin9_MinImputed.ARHGEF26(WT)vsMock.SMC.3xFlag.21JUL2020.tsv				
# Genoppi_Whitehead_martin9_MinImputed.HDAC9vsMock.EC.3xFlag.21JUL2020.tsv				
# Genoppi_Whitehead_martinBroad_MinImputed.ADAMTS7vsMock.SMC.expanded.21JUL2020.tsv				
# Genoppi_Whitehead_martinBroad_MinImputed.ARHGEF26vsMock.EC.nonSGS.21JUL2020.tsv				
# Genoppi_Whitehead_martinBroad_MinImputed.EDN1vsMock.SMC.expanded.21JUL2020.tsv				

## we expect two more rows (20 rows)
stopifnot(sum(summary$ok) == 20)
summary <- summary[summary$ok, ]
summary <- summary[order(summary$Cell.line),]
mypaths <- summary$Data.path
mypaths <- gsub('run055', 'run056', mypaths)
write.table(mypaths, file = 'run056_filtered_tier1_paths.tsv', sep = '\t')

## write table with IPs 
library(xlsx)
outfile = paste0(run,'_supplementary_table_1.xlsx')
options(java.parameters = "-Xmx8000m") # otherwise java is out of memory

## write to xlsx using java
for (i in seq_along(mypaths)){
  
  gc()
  
  # get file and naming convention
  p = mypaths[i]
  df = read.csv(p, sep = '\t')
  id = as.data.frame(t(as.matrix(unlist(strsplit(basename(p), '\\.|vs'))[1:5])))
  colnames(id) <- c('facility', 'submission', 'bait', 'control', 'cell')
  rownames(id) <- NULL
  sheetname = paste0(id$submission,'.',id$cell,'.',id$bait)
  write(sheetname, stdout())
  xlsx::write.xlsx(df, file = outfile, sheetName = sheetname, append = i > 1, row.names = F)

}




