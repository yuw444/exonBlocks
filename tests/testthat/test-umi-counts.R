library(data.table)

temp <- data.table::fread(
    "/scratch/g/chlin/Yu/exonBlocks/sub/raw/PBMCs.allruns.umicounts_intronexon.txt",
    header = TRUE,
    nrows = 1000
)
