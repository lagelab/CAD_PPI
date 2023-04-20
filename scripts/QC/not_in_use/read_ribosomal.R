
library(data.table)
library(readxl)

path <- "~/Downloads/Supplementary Table 4 - IP-MS results.xlsx"

# read all sheets of path into a list
all_sheets <- readxl::excel_sheets(path)
lst <- lapply(all_sheets, function(sheet) read_excel(path, sheet = sheet))
names(lst) <- all_sheets

# get ribosomal proteins for each and count how many are significant
counts <- do.call(rbind, lapply(lst, function(df){
    df <- df[grepl(df$gene, pattern = "(RPL)|(RPS)"),]
    data.frame(count=nrow(df), significant=nrow(df[df$significant == TRUE,]))
}))

counts$id <- rownames(counts)

# counts from supplementary table 5 (so that order is correct)
theorder <- scan(what=character())

# write and insert 
write.table(counts[match(theorder, counts$id),],    file="matched_counts.txt", sep="\t")



