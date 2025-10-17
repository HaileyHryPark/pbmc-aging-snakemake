library(devtools)
library(withr)

if(!file.exists(snakemake@output[[1]]))
  with_libpaths(new = snakemake@input, install_github(snakemake@params[["path"]]))
