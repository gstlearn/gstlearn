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
