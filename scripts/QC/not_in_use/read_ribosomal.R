
library(data.table)
library(readxl)

path <- "~/Downloads/Supplementary Table 4 - IP-MS results.xlsx"

# read all sheets of path into a list
all_sheets <- readxl::excel_sheets(path)
lst <- lapply(all_sheets, function(sheet) read_excel(path, sheet = sheet))
names(lst) <- all_sheets

# get ribosomal proteins for each and count how many are significant
counts <- do.call(rbind, lapply(lst, function(df){
    n <- nrow(df)
    n_all_sig <- nrow(df[df$significant == TRUE, ])
    df <- df[grepl(df$gene, pattern = "(RPL)|(RPS)"), ]
    n_ribosomal <- nrow(df)
    n_ribosomal_sig <- nrow(df[df$significant == TRUE, ])
    pct_ribosomal_sig <- n_ribosomal_sig/n_all_sig
    data.frame(count = n_ribosomal, significant = n_ribosomal_sig, pct_ribosomal_sig)
}))

counts$id <- rownames(counts)

# create a vector of the following strings:
the_order <- c("mB.EC.ARHGEF26", "m9.EC.HDAC9", "m11.EC.FN1", "m12.EC.EDN1",
    "m13.EC.PHACTR1", "m15.EC.FLT1", "m15.EC.PLPP3", "m17.EC.BCAS3",
    "m17.EC.JCAD", "mB.SMC.ADAMTS7", "mB.SMC.EDN1", "m5.SMC.BCAS3",
    "m6.SMC.HDAC9", "m9.SMC.ARHGEF26", "m10.SMC.FN1", "m10.SMC.JCAD",
    "m12.SMC.PHACTR1", "m13.SMC.PLPP3", "m14.SMC.EDNRA", "m15.SMC.ARHGEF26")


counts

# write and insert
counts <- counts[match(theorder, counts$id),]
write.table(counts, file = "~/Downloads/matched_counts.txt", sep = "\t")


getwd()
