

library(data.table)
setwd('/Users/flassen/Projects/14_micom_clean/MICOM/')

# test expected dataset and origin dataset
d1 <-fread("data/ms_data_genoppi_input/run056/Broad.mB.ARHGEF26vsMock.EC.nonSGS.20NOV2020.tsv")
d2 <- fread("data/mass_spec_broad/ARHGEF26/ARHGEF26_EC_nonSGS.csv")
nrow(d1) # 1930
nrow(d2) # 1963

d1 <- d1[,c('gene','rep1','rep2')]
d2 <- d2[,c('gene','rep1', 'rep2')]

# some genes are duplicated and then when we
# match replicate FC we get a different answer,
# thus, let's remove them and match.
dups <- c(d1$gene[duplicated(d1$gene)],
          d2$gene[duplicated(d2$gene)])

d1 <- d1[!d1$gene %in% dups,]
d2 <- d2[!d2$gene %in% dups,]

# merging replciate FC
mrg <- merge(d1,d2, by = "gene")
nrow(mrg)
sum(mrg$rep1.x == mrg$rep1.y)/nrow(mrg) # 1
sum(mrg$rep2.x == mrg$rep2.y)/nrow(mrg) # 1

