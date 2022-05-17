# Add automatic display for all AStringable objects (see onAttach comment below)
setMethod(f = "show", signature = "_p_AStringable", definition = function(object){ AStringable_display(object) })

# Add [] set/get operator to Db calss
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
