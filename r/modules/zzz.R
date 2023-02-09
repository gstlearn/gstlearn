#
# Initialization of plot.R
#
".onAttach" <- 
function(libname, pkgname)
{
	pos = match(paste0("package:", pkgname), search())
	plot.initialize(pos)
}
