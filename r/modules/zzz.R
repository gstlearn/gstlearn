# Add automatic display for all AStringable objects (see onAttach comment below)
setMethod(f = "show", signature = "_p_AStringable", definition = function(object){ AStringable_display(object) })

# Add [] accessors to Db class
setMethod("["  ,signature(x="_p_Db"),
  function(x,i,j,...,drop) 
  { 
    if (! hasArg(j))
      x$getItem(i)
    else
      x$getItem(i,j)
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

".onAttach" <- 
function(...)
{
  # Cannot call setMethod in onAttach (do it outside !)
  
  packageStartupMessage("*** Welcome to the gstlearn package ***\n\n")
}
