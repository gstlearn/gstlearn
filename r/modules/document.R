################################################################################
#                                                                              #
#                            gstlearn R package                                #
#                                                                              #
# Copyright (c) (2023) MINES Paris / ARMINES                                   #
# Authors: gstlearn Team                                                       #
# Website: https://gstlearn.org                                                #
# License: GPL v3                                                              #
#                                                                              #
################################################################################

# This file contains some helping functions used specifically for documentation
# of gstlearn package for R language.
#

# The various pieces of documentation are supposed to be located
# at the following URL
urlMP = "https://soft.minesparis.psl.eu/gstlearn"

#' Check if Internet is available
#' This function requires the package 'lares' to be installed
#' @return TRUE if Internet is available and FALSE otherwise
isInternetAvailable <- function()
{
  flag = FALSE
  if (require("lares", quietly=TRUE))
    flag = haveInternet()
  flag
}

#' Display a piece of documentation (Markdown file) from 'references' directory
#' @param filename Name of the file containing the text to be displayed
displayMarkdown <- function(filename)
{
  if (isInternetAvailable())
    cat(readLines(paste0(c(urlMP, "references", filename), collapse='/'), warn=FALSE), sep="\n")
  else
    cat(readLines(filename, warn=FALSE), sep="\n")
  invisible()
}

#' Load and returns the contents of the data file 'filename' located within 'data/directory'
#' @param directory Name of the Directory located in the 'data' sub-directory of the URL
#' @param filename Name of the data file to be loaded
#' @return The name of the returned data file
loadData <- function(directory, filename)
{
  if (isInternetAvailable())
    download.file(paste0(c(urlMP, "data", directory, filename), collapse='/'), filename, quiet=TRUE)
  else
    filename = file.path('.', "data", directory, filename)
  filename
}