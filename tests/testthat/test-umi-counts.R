library(data.table)

temp <- data.table::fread(
    "/home/yu-wang/Documents/exonBlocks/meta/features/PBMCs.allruns.umicounts_intronexon.txt",
    header = TRUE,
    nrows = 1000
)
