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

# The various pieces of documentation are supposed to be located at the following URL
urlGST = "https://soft.minesparis.psl.eu/gstlearn"

# Next definitions are used to decorate the MD files displayed in rstudio
header = c(
  "<style>md-block { color:gray; background-color:white; }</style>",
  "<md-block>\n")
trailer = c(
  "</md-block>")
  
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

#'
#' Return the complete name of a file:
#' - if Internet is available, the file is retrieved from the web site
#' - if not, it is assumed to be present locally
#'    
#' @param directory: Name of the target directory (used for 'where' = "data", "None" otherwise)
#' @param filename: Name of the file to be downloaded
#' @param where: 'data' or 'graphics' or 'mdfile'
#'    
#' @remarks
#' When retrieving file, and according to 'where', the origin is:
#' - when 'where' == "graphics"
#'   urlGST + "/references' + '/Figures' + '/'
#' - when 'where' == "mdfile"
#'   urlGST + "/references' + '/'
#' - when 'where' == "data"
#'   urlGST + "/data" + '/' + directory + '/' + filename
downloadRemoteFile <- function(directory, filename, where)
{
  # Generate local name

  if (where == 'graphics')
  {
    localname = paste0(c('Figures' ,filename), collapse='/')
    if (!dir.exists('Figures')) {dir.create('Figures')}
  }
  else if (where == 'mdfile')
    localname = paste0(c('.' ,filename), collapse='/')
  else if (where == 'data')
    localname = paste0(c('.' ,filename), collapse='/')
  else
  {
    print("'downloadRemoteFile' does not know about 'where' = ", where)
    exit()
  }

  if (isInternetAvailable())
  {
    if (where == 'graphics')
        pathname = paste0(c(urlGST, 'references',  'Figures', filename), collapse='/')
    else if (where == 'mdfile')
        pathname = paste0(c(urlGST, 'references', filename), collapse='/')
    else if (where == 'data')
        pathname = paste0(c(urlGST, "data", directory, filename), collapse='/')
    else
        print("'downloadRemoteFile' does not know about 'where' = ", where)
    
    # The file is loaded in the local environment
    err = download.file(pathname, localname, quiet=TRUE)
  }
  localname
}

loadFigure <- function(filename)
{
  downloadRemoteFile(None, filename, "graphics")
}

#' Returns the decorated documentation (Markdown file) from 'references' directory
#' @param filename Name of the file containing the text to be displayed
#' TODO: the color does not function... to be fixed.
#' remark: the returned string must be displayed in a RMarkdown chunk as follows:
#'   {r, echo=FALSE, result='asis'}
#'    cat(XXX, sep="\n")
#'   }
loadDoc <- function(filename)
{
  filepath = downloadRemoteFile(NULL, filename, "mdfile")
  multiline = readLines(filepath, warn=FALSE)
  
  # Loop on the lines to detect graphic
  searchItem = "(Figures/"
  for (i in 1:length(multiline))
  {
    targetLine = multiline[i]
    if (grepl(searchItem, targetLine, fixed=TRUE))
    {
      graphicFile = sub(".*Figures/", "", targetLine)[1]         # Extract file name
      graphicFile = substr(graphicFile, 1, nchar(graphicFile)-1) # suppress last character
      filefig = downloadRemoteFile(NULL, graphicFile, "graphics")
    }
   }
  result = c(header, multiline, trailer)
  result
}

#' Load and returns the contents of the data file 'filename' located within 'data/directory'
#' @param directory Name of the Directory located in the 'data' sub-directory of the URL
#' @param filename Name of the data file to be loaded
#' @return The name of the returned data file
loadData <- function(directory, filename)
{
  downloadRemoteFile(directory, filename, "data")
}

#'
#' This function is used to clean the files loaded for find_statement_documentation
#' @param filename: Name of the target file
cleanDoc <- function(filename)
{
  # Remove the target file
  if (file.exists(filename)) 
  file.remove(filename)

  # Remove the downloaded graphic files (in subdirectory 'Figures')
  if (dir.exists('Figures'))
    unlink('Figures', recursive=TRUE)
}
