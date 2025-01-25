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
#
# Initialization of plot.R
#
".onAttach" <- 
function(libname, pkgname)
{
	#OptDbg_reset() # Remove this as it is obvious and it makes roxygen crashing
	plot.initialize()
}

".onLoad" <- function(libname, pkgname)
{

getMethodsOfClass = function(name)
{
  x=ls(envir=asNamespace("gstlearn"))
  split_vec = strsplit(x,"_")
  first_part <- sapply(split_vec, function(x) x[1])
  second_part <- sapply(split_vec, function(x) x[2])
  len = sapply(split_vec,function(x) length(x))
  ind = (first_part == name) & len == 2
  return(second_part[ind])
}

addMethods = function(base, derived) {
  setMethod(
    "$", paste0("_p_", derived),
    function(x, name) {
      accessorFuns = list()
      for (classe in c(base, derived))
      {
        methodsName = getMethodsOfClass(classe)
        for (namemethod in methodsName)
        {
          accessorFuns[[namemethod]] <- get(paste0(classe, "_", namemethod))
        }
      }
      idx = match(name, names(accessorFuns))
      if (is.na(idx)) {
        return(callNextMethod(x, name))
      }
      f = accessorFuns[[idx]]
      result = function(...) {
        f(x, ...)
      }
      return(result)
    }
  )
}

addMethods("ModelCovList", "ModelGeneric")
}
