library("knitr")
library("callr")

# Retrieve command args
args   <- commandArgs(trailingOnly = TRUE)
script <- args[1]
outdir <- args[2]

scriptname = basename(script)
scriptname = tools::file_path_sans_ext(scriptname)

outscript = file.path(outdir, paste0(scriptname, ".R"))
outpath   = file.path(outdir, paste0(scriptname, ".out"))

# 1. Convert Rmd to R script (no documentation)
# https://stackoverflow.com/questions/71183578/how-to-extract-all-code-from-an-rmarkdown-rmd-file
knitr::purl(input = script, output = outscript, documentation = 0)

# 2. Execute R script and dump output in a text file
# https://callr.r-lib.org/reference/rscript.html
rscript(outscript, stdout=outpath)


