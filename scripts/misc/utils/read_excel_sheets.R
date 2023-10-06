

# load data for enrichment
read_excel_sheets <- function(path){
  sheets <- readxl::excel_sheets(path)
  data <- lapply(sheets, function(sheet) readxl::read_excel(path, sheet))
  names(data) <- sheets
  return(data)
}