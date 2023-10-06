setwd('~/Projects/03_MICOM/')

# ready paths and prefixes
date = '201120'
run = 'run056'
target = 'derived/maps/maps-run056/'

# setup paths and run data
paths = read.csv('run056_filtered_tier1_paths.tsv', sep = '\t')$x
data = list()

# what baits and what cells?
baits = unique(unlist(lapply(paths, function(x) unlist(strsplit(x, split = '\\.|vs'))[3])))
cells = unique(unlist(lapply(paths, function(x) unlist(strsplit(x, split = '\\.|vs'))[5])))
re_cells = c(cells, paste(cells, collapse = '|'))

for (cell in re_cells){
  for (bait in baits){
    
    # get interactors by bait and cell
    res = paths[grepl(cell, paths) & grepl(bait, paths)]
    res_data = lapply(res, function(x) read.csv(x, sep = '\t'))
    res_data = do.call(rbind, res_data)
    
    if (!is.null(res_data)){
      
      # filter out recurrent interactors or ambigious ones
      res_data = res_data[,c(1,7,11,12,13)]
      res_data = res_data[!duplicated(res_data),]
      ambigious_interactor = duplicated(res_data$gene) & !res_data$significant
      res_data = res_data[!ambigious_interactor,]
      
      # append conditions to data (Yu-Hans convention)
      cell_string = gsub('EC\\|SMC','EC-SMC',cell)
      res_data$ListName <- paste0(bait,'_',cell_string)
      res_data$Bait <- bait
      res_data$CellType <- cell_string
      
      # what file to write
      outfile = file.path(target, paste0(run,'_',bait,'_',cell_string,'.InteractorTable.txt'))
      write(outfile, stdout())
      
      # reorder columns and write file
      res_data = res_data[,c(6:8,1:5)]
      colnames(res_data) <- c('ListName',	'Bait',	'CellType',	'GeneName',	'Interactor',	'Chr',	'Start',	'End')
      write.table(res_data, outfile, sep = '\t', quote = F, row.names = F)
      data[[outfile]] <- res_data
      
    }
  }
}

# get combined tables and write
#table_ec <- do.call(rbind, data[grepl('EC\\.',names(data))])
#table_smc <- do.call(rbind, data[grepl('SMC\\.',names(data))])
#table_ecsmc <- do.call(rbind, data[grepl('EC-SMC\\.',names(data))])
#prefix <- paste0('derived/maps/',date)
#write.table(table_ec, file = prefix(), sep = '\t')

## end