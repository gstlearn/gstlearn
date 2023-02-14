library("rmarkdown")
library("knitr")

# Retrieve command args
args   <- commandArgs(trailingOnly = TRUE)
script <- args[1]
outdir <- args[2]

outfile = basename(script)
outfile = tools::file_path_sans_ext(outfile)
outfile = paste0(outfile, ".md")

outpath = file.path(outdir, outfile)

#cat("Running  script:", script, "\n")
#cat("         outdir:", outdir, "\n")
#cat("        outpath:", outpath, "\n")

# Execute Rmd in a new environment
# https://stackoverflow.com/questions/53695596/r-markdown-markdown-workspace-to-r-workspace
rmarkdown::render(script, quiet=TRUE, envir = new.env(),
		          output_format="md_document", output_file=outfile, output_dir=outdir)

# Remove image file absolute paths from the output file = all lines starting with 
#-     ![](     i.e. :
# ![](/home/fors/Projets/gstlearn/gstlearn/build/tests/rmd/Release/Starting_files/figure-markdown_strict/unnamed-chunk-3-1.png)
#-     <img     i.e. :
# <img src=/home/fors/Projets/gstlearn/gstlearn/build/tests/rmd/Release/Starting_files/figure-markdown_strict/unnamed-chunk-3-1.png>
# https://stackoverflow.com/questions/41575419/text-mining-in-r-remove-rows-from-text-file-starting-with-keywords

#cat("Patching output file:", outpath, "\n")

df = readLines(outpath)
x = grepl("^!\\[\\]\\(", df)
df <- df[!x]
x = grepl("^<img", df)
df <- df[!x]
writeLines(df, outpath)
