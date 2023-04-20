# get interactor statistics
# by frhl
# data: 21-04-19

library(readxl)

# what genesets are enriched in the WT versus MT.
setwd('~/Projects/14_micom_clean/MICOM/')

path <- 'derived/MICOM_DataSummary.xlsx'
readxl::excel_sheets(path)

V <- list(
  EC = read_xlsx(path, )
  SMC = fread('data/interactor_lists/run056/COMBINED_SMC.InteractorTable.txt'),
  EC_SMC = fread('data/interactor_lists/run056/COMBINED_EC-SMC.InteractorTable.txt')
)

lapply(V, function(x){
  
})
