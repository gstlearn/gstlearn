# Example of libdir :
# /home/fors/Projets/gstlearn/gstlearn/build/r/Release/gstlearn

load_gstlearn = function(libdir)
{
  dyn.load(paste0(libdir,"/gstlearn.so"))
  cacheMetaData(1)
  gstlearn <- new.env()
  eval(parse(paste0(libdir,"/gstlearn.R")), envir=gstlearn)
  attach(gstlearn)

  # Add default show S4 methods to all objects (assuming that they inherits from AStringable)
  cl = showMethods('$', printTo=FALSE)
  cl = cl[grep("^x=\"_p_", cl)]
  cl = gsub("^x=\"","",cl)
  cl = gsub("\"","",cl)
  invisible(lapply(cl,function(x) { 
    cmd = paste0("setMethod(f = \"show\", signature = \"",x,"\"",
                 ", definition = function(object){ ",
                 "AStringable_display(object) })")
    eval(parse(text=cmd))
  }))
  
  
  source(paste0(libdir,"/plot.r"))
  
  setMethod("["  ,signature(x="_p_Db"),
  	function(x,i,j,...,drop) 
  	{ 
  		if (! hasArg(j))
  		{
	  		ans = x$getItem(i)
	  		names = x$getItemNames(i)
	  	}
	  	else
	  	{
	  		ans = x$getItem(i,j)
	  		names = x$getItemNames(j)
	  	}
	  	
	  	# Conversion to numeric vector (dim=1) or data.frame (dim>1)
	  	# is performed here. It would be beneficial to perform
	  	# this conversion in swig R-dependent module	
	  	if (length(ans) == 1) 
	  		ans = as.numeric(ans[[1]])
	  	else
	  	{
	  		ans = as.data.frame(ans)
	  		names(ans) <- names
	  	}
 		ans
  	}
  )
  setMethod("[<-",signature(x="_p_Db"),
  	function(x,i,j,value)
  	{
  		if (! hasArg(j))
  			x$setItem(i,value)
  		else
  			x$setItem(i,j,value)
  		x
  	}
  )
}
