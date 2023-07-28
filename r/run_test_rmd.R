library("knitr")
library("callr")

# Retrieve command args
args   <- commandArgs(trailingOnly = TRUE)
script <- args[1]
outdir <- args[2]
out_type = "R"
if (length(args) > 2) {
  out_type = args[3]
}

scriptname = basename(script)
scriptname = tools::file_path_sans_ext(scriptname)

outscript = file.path(outdir, paste0(scriptname, ".", out_type))
outpath   = file.path(outdir, paste0(scriptname, ".out"))

cat("Sys.getenv('GSTLEARN_DATA')=",Sys.getenv('GSTLEARN_DATA'),"\n")

if (out_type == 'R') {
  # 1. Convert Rmd to R script (no documentation)
  # https://stackoverflow.com/questions/71183578/how-to-extract-all-code-from-an-rmarkdown-rmd-file
  knitr::purl(input = script, output = outscript, documentation = 0)

  # 2. Execute R script and dump output in a text file
  # https://callr.r-lib.org/reference/rscript.html
  rscript(outscript, stdout=outpath)
} else if (out_type == "html") {
  rmarkdown::render(script, output_format='html_document', quiet=FALSE, output_dir=outdir)
} else {
  stop("Hun ?\n")
}

