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
	pos = match(paste0("package:", pkgname), search())
	OptDbg_reset()
	plot.initialize(pos)
}
