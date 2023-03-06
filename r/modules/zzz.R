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
