library(readxl)



reactome <- list(
  reactome = 'derived/210422_reactome_network_only_table.xlsx',
  hallmark = 'derived/210422_msigdb_h_network_only_table.xlsx'
)

go <- list(
  mf = 'derived/210422_mf_network_only_table.xlsx',
  bp = 'derived/210422_bp_network_only_table.xlsx',
  cc = 'derived/210422_cc_network_only_table.xlsx'
)




combine <- function(lst){
  
  innames <- names(lst)
  result <- lapply(innames, function(cur_name){
    cur_path <- lst[[cur_name]]
    cur_df <- read_excel_sheets(cur_path)
    names(cur_df) <- paste0(cur_name,'.',names(cur_df))
    return(cur_df)
  })
  return(unlist(result, recursive = F))
}




write_xlsx(combine(reactome), path = 'Supplementary Table 10 - 210422_reactome_only_table.xlsx')
write_xlsx(combine(go), path = 'Supplementary Table 111 - 210422_GO_only_table.xlsx')



