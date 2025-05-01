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
package_version = packageVersion("gstlearn")

# Next definitions are used to decorate the MD files displayed in rstudio
header = c(
  "<style>md-block { color:gray; background-color:white; }</style>",
  "<md-block>\n")
trailer = c(
  "</md-block>")
  
#' Check if Internet is available
#' This function requires the package 'lares' to be installed
#' @return TRUE if Internet is available and FALSE otherwise
internetAvailable <- function()
{
  flag = FALSE
  if (require("lares", quietly=TRUE))
    flag = haveInternet()
  flag
}

#' Return the absolute path of a file:
#' - it is assumed to be present locally in '.' ('where' and 'directory' are ignored)
#' - if not, it is assumed to be present locally in './doc/<where>', '../../doc/<where>' or '../../<where>'
#' - if not, if the GSTLEARN_DIR environment variable is defined, it is assumed to be present in '<GSTLEARN_DIR>/gstlearn/doc/<where>'
#' - if not, if Internet is available, the file is downloaded from the gstlearn website in a temporary file
#' 
#' @param filename: Name of the file to be located
#' @param where: 'data' or 'references'
#' @param directory: Name of the data file directory (only used for 'where' = "data")
#' @param verbose: TRUE to activate verbose mode
#' @param version: Use a specific gstlearn version when searching the file on the web (string)
#' 
locateFile <- function (filename, where='references', directory=NULL, verbose=FALSE, version=package_version)
{
  argfilename = filename
  if (verbose)
    print(paste("Current directory is", getwd()))
  
  # Test current directory
  localname = file.path('.', filename)
  if (file.exists(localname))
  {
    fullname = normalizePath(localname)
    if (verbose)
      print(paste(filename, "found... Full path is", fullname))
    return(fullname)
  }
  else if (verbose)
  {
    print(paste(localname, "not found..."))
  }
  
  # Test locally in other directories
  if (!(where %in% c('references', 'data')))
  {
    print(paste("'locateFile' does not know about 'where' = ", where))
    return(NULL)
  }
  if (where == 'data' && !is.null(directory))
    filename = file.path(directory, filename)
  
  folders = list(file.path('.',"doc",where),
                 file.path('..','..',"doc",where),
                 file.path('..','..',where))
  for (f in folders)
  {
    localname = file.path(f, filename)
    if (file.exists(localname))
    {
      fullname = normalizePath(localname)
      if (verbose)
        print(paste(filename, "found... Full path is", fullname))
      return(fullname)
    }
    else if (verbose)
    {
      print(paste(localname, "not found..."))
    }
  }

  # Test in GSTLEARN_DIR environment variable
  if (!is.null(Sys.getenv("GSTLEARN_DIR")) &&
               Sys.getenv("GSTLEARN_DIR") != "")
  {
    localname = file.path(Sys.getenv("GSTLEARN_DIR"), "gstlearn", "doc", where, filename)
    if (file.exists(localname))
    {
      fullname = normalizePath(localname)
      if (verbose)
        print(paste(filename, "found... Full path is", fullname))
      return(fullname)
    }
    else if (verbose)
    {
      print(paste(localname, "not found..."))
    }
  }

  # Test on the web
  if (!internetAvailable())
  {
    print(paste("Error: Cannot access to", filename, "(no Internet)!"))
    return(NULL)
  }
  
  # Download from Internet in a temporary file
  localname = paste0(urlGST, '/', version, '/', where, '/', directory, '/', argfilename)
  fullname = tempfile()
  if (!download.file(localname, fullname, quiet=TRUE))
  {
    if (verbose)
      print(paste(localname, "found... Full path is", fullname))
    return(fullname)
  }

  print(paste("Cannot access URL:", localname, "!"))
  return(NULL)
}

#' Returns the decorated documentation (Markdown file) from 'references' directory
#' 
#' @param filename Name of the Markdown file containing the text to be displayed
#' @param verbose: TRUE to activate verbose mode
#' @param version: Use a specific gstlearn version when searching the file on the web (string)
#' 
#' TODO: the color does not function... to be fixed.
#' remark: the returned string must be displayed in a RMarkdown chunk as follows:
#'   {r, echo=FALSE, result='asis'}
#'    cat(XXX, sep='new_line')
loadDoc <- function(filename, verbose=FALSE, version=package_version)
{
  if (!require("stringr", quietly=TRUE))
  {
    print(paste("'stringr' package must be installed!"))
    return("")
  }
    
  filepath = locateFile(filename, verbose=verbose, version=version)
  if (is.null(filepath))
    return(paste("File ", filename, "not found!"))
  
  multiline = readLines(filepath, warn=FALSE)
  
  # Capture Markdown images (beginning ![description](filename) ending)
  pattern = '(.*)\\!\\[(.*)\\]\\((.+)\\)(.*)'
  for (i in 1:length(multiline))
  {
    targetLine = multiline[i]
    img = str_match(targetLine, pattern)
    if (!is.na(img[1]))
    {
      beginning = img[2]
      imgdesc = img[3]
      imgfile = locateFile(img[4], verbose=verbose, version=version)
      ending = img[5]
      if (is.null(imgfile))
        return(paste("File", img[4], "not found!"))
      # Convert in base64 for embedding the image (not needed n R)
      #with open(imgfile, 'rb') as image_file: # this is python code
      #    imgfile = base64.b64encode(image_file.read())
      # Reconstruct the full Markdown line
      multiline[i] = paste0(beginning,'![',imgdesc,'](',imgfile,')',ending)
    }
  }
  result = c(header, multiline, trailer)
  result
}

#' Returns path of a data file 'filename' located within the 'data' directory (locally or from the web site)
#' @param directory Name of the sub-irectory located in the 'data' directory
#' @param filename Name of the data file to be loaded
#' @param verbose: TRUE to activate verbose mode
#' @param version: Use a specific gstlearn version when searching the file on the web (string)
#' @return The name of the returned data file
loadData <- function(directory, filename, verbose=FALSE, version=package_version)
{
  locateFile(filename, "data", directory, verbose=verbose, version=version)
}
