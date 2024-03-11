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

if (out_type == 'R') {
  
  # Load dependency
  library("callr")
  
  # 1. Convert Rmd to R script (no documentation)
  # https://stackoverflow.com/questions/71183578/how-to-extract-all-code-from-an-rmarkdown-rmd-file
  knitr::purl(input = script, output = outscript, documentation = 0)

  # 2. Execute R script and dump output in a text file
  # https://callr.r-lib.org/reference/rscript.html
  rscript(outscript, stdout=outpath)
  
  # 3. Remove unwanted logs
  tx  <- readLines(outpath)
  # Remove image paths included by MD files
  # ![Description](/home/drenard/project_gstlearn/gstlearn/build/doc/references/Figures/uniformdiscretpmf.png
  tx2  <- gsub(pattern = ".*/Figures/.*", replace = "XXX", x = tx)
  writeLines(tx2, con=outpath)
  
} else if (out_type == "html") {
  
  # Load dependency
  library("knitr")
  
  # Execute Rmd script and dump output in a html file (self_contained see html header in Rmd)
  rmarkdown::render(script, output_format='html_document', quiet=FALSE, output_dir=outdir)
  
} else {
  
  stop("Hun ?\n")
  
}

