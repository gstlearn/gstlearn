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
# This is a set of functions working with ggplot2() which enable performing
# plots easily.
# It also requires the packages:
# - ggnewscale: used to reset color scale
# - ggrepel: used to take care of placing the labels to minimize overlap 
#

#' Define the global values in the given environment position (search)
#
#' @param pos Position in the list of packages
plot.initialize <- function(pos=1) 
{
  assign("plot.defaultDims", list(c(8,8), c(8,8)), pos=pos)
  assign("plot.defaultXlim", list(c(NA,NA), c(NA,NA)), pos=pos)
  assign("plot.defaultYlim", list(c(NA,NA), c(NA,NA)), pos=pos)
  assign("plot.defaultAspect",  c(0, 1), pos=pos)
  invisible()
}

#' Check if an argument is defined
#' @param arg Argument to be checked
#' @noRd
.isNotDef <- function(arg)
{
  if (is.null(arg)) return (TRUE)
  warn.old = options("warn")
  options(warn = -1)
  if (length(arg) == 1)
  {
    if (is.na(arg)) 
    {
      options(warn.old)
      return (TRUE)
    }
  }
  else
  {
    for (i in 1:length(arg))
    {
      if (is.na(arg[i])) 
      {
        options(warn.old)
        return (TRUE)
      }
    }
  }
  options(warn.old)
  return (FALSE)
}

#' Set the default values for all subsequent Geographical figures
#'
#' @param dims Vector giving the dimensions of the figure
#' @param xlim Bounds of the figure along the horizontal axis (when left to NA, it will be adjusted to the figure contents)
#' @param ylim Bounds of the figure along the vertical axis (when left to NA, it will be adjusted to the figure contents)
#' @param asp Aspect ratio Y/X
plot.setDefaultGeographic <- function(dims=NA, xlim=NA, ylim=NA, asp=NA)
{
  .plot.setDefaultInternal(2, dims=dims, xlim=xlim, ylim=ylim, asp=asp)
}

#' Set the default values for all subsequent non-Geographical figures
#'
#' @param dims Vector giving the dimensions of the figure
#' @param xlim Bounds of the figure along the horizontal axis (when left to NA, it will be adjusted to the figure contents)
#' @param ylim Bounds of the figure along the vertical axis (when left to NA, it will be adjusted to the figure contents)
#' @param asp Aspect ratio Y/X
plot.setDefault <- function(dims=NA, xlim=NA, ylim=NA, asp=NA)
{
  .plot.setDefaultInternal(1, dims=dims, xlim=xlim, ylim=ylim, asp=asp)
}

#' Set the default values for all subsequent Geographical figures
#' @param mode 1 for Geographical and 2 for non-Geographical parameteres
#' @param dims Dimensions of the figures
#' @param xlim Bounds along the horizontal axis
#' @param ylim Bounds along the vertical axis
#' @param asp Aspect ratio Y/X
#' @noRd
.plot.setDefaultInternal <- function(mode=1, dims=NA, xlim=NA, ylim=NA, asp=NA)
{
  if (!.isNotDef(dims))
    plot.defaultDims[[mode]] = dims
  if (!.isNotDef(xlim))
    plot.defaultXlim[[mode]] = xlim
  if (!.isNotDef(ylim))
    plot.defaultYlim[[mode]] = ylim    
  if (!.isNotDef(asp))
    plot.defaultAspect[[mode]] = asp
}

#' Print the Default values for both Geographical and non-geographical subsequent figures
plot.printDefault <- function()
{
  for (mode in 1:2)
  {
    if (mode == 1)
      cat("Non geographical defaults (mode=1):\n")
    else
      cat("Geographical defaults (mode=2):\n")
    
    if (!.isNotDef(plot.defaultDims[[mode]]))
      cat("- Figure dimensions =", plot.defaultDims[[mode]],"\n")
    else
      cat("- Figure dimensions (not defined)\n")
    
    if (!.isNotDef(plot.defaultXlim[[mode]]))
      cat("- Limits along X =",plot.defaultXlim[[mode]],"\n")
    else
      cat("- Limits along X (not defined)\n")
    
    if (!.isNotDef(plot.defaultYlim[[mode]]))
      cat("- Limits along Y =",plot.defaultYlim[[mode]],"\n")
    else
      cat("- Limits along Y (not defined)\n")
    
    if (plot.defaultAspect[mode] != 0)
      cat("- Aspect =",plot.defaultAspect[mode],"\n")
    else
      cat("- Aspect (automatic)\n")
  }  
}

#' Define the color map
#' @param palette Reference palette used for defining the current color map
#' @noRd
.scaleColorFill <- function(palette, naColor=NA, ...)
{
  rcb <- c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", "PuBu",
      "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd",
      "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3",
      "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")
  v <- c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo",
      "A", "B", "C", "D", "E", "F", "G", "H")
  rcb_num <- 1:18
  
  aes_list = c("colour", "fill")
  if (length(palette) == 0) 
  {
    layer = scale_color_gradient(na.value=naColor, ...)
  }
  else if(length(palette) == 1) 
  {
    if (any(palette == rcb) | any(palette == rcb_num)) 
    {
      layer = scale_color_distiller(palette=palette, aesthetics=aes_list, 
 	     na.value=naColor, ...)
    } 
    else if(any(palette == v)) 
    {
      layer = scale_color_viridis_c(option=palette, aesthetics=aes_list, 
 	     na.value = naColor, ...)
    } 
  } 
  else if(length(palette) == 2) 
  {
    low = palette[1]
    high = palette[2]
    layer = scale_color_gradient(low= low, high= high, aesthetics=aes_list, 
 	   na.value=naColor, ...)
  } 
  else if(length(palette) == 3) 
  {
    low = palette[1]
    mid = palette[2]
    high = palette[3]
    layer = scale_color_gradient2(low= low, mid= mid, high= high, 
 	   aesthetics=aes_list, na.value=naColor, ...)
  } 
  else 
  {
    layer = scale_colour_manual(values= palette, aesthetics=aes_list, 
 	   na.value = naColor, ...)
  }
  layer
}

#' Define a series of distinct colors
#' @noRd
.getColors <- function()
{
  c("blue", "red", "green", "brown", "orange", "purple", "yellow")
}

#' Print the contents of a ggplot, possibly without warnings
#'
#' @param p Current contents of the ggplot()
#' @param flagSuppressWarnings TRUE to suppress informational warnings
ggPrint <- function(p, flagSuppressWarnings = TRUE)
{
  if (flagSuppressWarnings)
    suppressWarnings(plot(p))
  else
    plot(p)
  invisible()
}

#' Initiate a Geographical display (using the default values for parameters)
#' @param figsize Optional parameter giving the dimensions of the figure
#' @return The ggplot object
#' @note When 'figsize' is defined, it overwrites the default dimensions 
#' @note coming from the geographical and non-geographical environments
#' @note Use printDefault() to visualize them and setDefaultGeographic() to modify them
ggDefaultGeographic <- function(figsize=NA)
{
  p <- ggplot()
  mode = 2
  
  if (.isNotDef(figsize))
 	locdims = plot.defaultDims[[mode]]
  else
  	locdims = figsize
   
  p <- p + plot.geometry(dims=locdims, 
                         xlim=plot.defaultXlim[[mode]], 
                         ylim=plot.defaultYlim[[mode]], 
                         asp=plot.defaultAspect[mode])
  p
}

#' Initiate a non-geographical display (using the default values for parameters)
#' @return The initiated ggplot object
#' @param figsize Optional parameter giving the dimensions of the figure
#' @note Use printDefault() to visualize them and setDefault() to modify them
ggDefault <- function(figsize=NA)
{
  p <- ggplot()
  mode = 1
  
  if (.isNotDef(figsize))
 	locdims = plot.defaultDims[[mode]]
  else
  	locdims = figsize
  
  p <- p + plot.geometry(dims=locdims, 
                         xlim=plot.defaultXlim[[mode]], 
                         ylim=plot.defaultYlim[[mode]], 
                         asp=plot.defaultAspect[mode])
  p
}

#' Check if the argument can be considered as an array (with possibly required dimensions)
#' @param arg Input argument
#' @param ndim Required dimension for the input argument (no check is performed if NA)
#' @noRd
.isArray <- function(arg, ndim=NA)
{
  if (length(arg) <= 1) return (FALSE)
  if (.isNotDef(arg)) return (FALSE)
  if (length(arg) != ndim) return (FALSE)
  
  TRUE
}

#' Draw the decoration of a figure (title, axis labels, ...)
#'
#' @param xlab Label along the horizontal axis
#' @param ylab Label along the vertical axis
#' @param title Title of the figure
#' @return The ggplot object
plot.decoration <- function(xlab = NA, ylab = NA, title = NA)
{
  p = list()
  if (!.isNotDef(xlab))
    p <- append(p, list(labs(x = xlab)))
  if (!.isNotDef(ylab))
    p <- append(p, list(labs(y = ylab)))
  if (!.isNotDef(title))
  {
    p <- append(p, list(labs(title = title)))
    p <- append(p, list(theme(plot.title = element_text(hjust = 0.5))))
  }
  p
}

#' Set the Geometry for the current plot
#'
#' @param dims Dimension of the figure
#' @param xlim Bounds along the horizontal axis
#' @param ylim Bounds along the vertical axis
#' @param asp  Aspect Ratio Y/X ("0" for an automatic aspect ratio)
#' @param expand Adding padding around data
#' @return The ggplot object
plot.geometry <- function(dims=NA, xlim=NA, ylim=NA, asp=NA, expand=waiver())
{
  p = list()
  if (! .isNotDef(dims[1]))
  {
    if (.isArray(dims, 2))
    {
      options(repr.p.width  = dims[1], repr.p.height = dims[2])
    }
    else
      cat("'dims' should be [a,b]. Ignored\n")
  }
  
  if (.isArray(xlim, 2))
  {
    p <- append(p, scale_x_continuous(limits=xlim, expand=expand))
  }

  if (.isArray(ylim, 2))
  {
    p <- append(p, scale_y_continuous(limits=ylim, expand=expand))
  }
  
  if (!.isNotDef(asp))
  {     
    if (asp != 0)
      p = append(p, coord_fixed(asp))
    else
      p = append(p, list(theme(aspect.ratio = 1)))
  }
  p
}

#' Function for representing a Model
#' @param model An object of class Model from gstlearn
#' @param ivar Rank of the variable to be represented (-1 for all variables)
#' @param jvar Rank of the second variable to be represented in multivariate case (-1 for all variables
#' @param vario An object of class Vario of gstlearn (optional)
#' @param idir Rank of the direction
#' @param ... Arguments passed to plot.varmod()
#'
#' If 'vario' is defined, the calculation direction is given by the definition of direction 'idir' within 'vario'
#' 
#' @return The ggplot object
plot.model <- function(model, ivar=0, jvar=0, vario=NA, idir=0, ...)
{
  p = list()
  p = append(p, plot.varmod(vario=vario, model=model, ivar=ivar, jvar=jvar, idir=idir,
  			 drawVario=FALSE, ...))
  p = append(p, plot.decoration(xlab = "Distance", ylab = "Variogram"))
  p  
}

#' Function for representing the Experimental Variogram
#' @param vario An object of the class Vario of gstlearn
#' @param ivar Rank of the variable to be represented 
#' @param jvar Rank of the second variable to be represented in multivariate case 
#' @param idir Rank of the direction to be represented 
#' @param ... Arguments passed to plot.varmod()
#' @return The ggplot object
plot.vario <- function(vario, ivar=-1, jvar=-1, idir=-1,...)
{
  p = list()
  p = append(p, plot.varmod(vario=vario, ivar=ivar, jvar=jvar, idir=idir,...))
  p
}

#' Select a list of items
#' @param nvalues Number of items in the list
#' @param sitem   Default value for the item withon this list (if non-negative)
#' @return The returned list of items 
.selectItemsInList <- function(nvalues, sitem=-1)
{
  if (sitem >= 0)
    outs = sitem
  else
    outs = seq(0, nvalues-1)
  outs
}

#' Draw an elementary experimental variogram
#' @noRd
.varioElementary <- function(vario, ivar=0, jvar=0, idir=0, 
    varColor='black', varLinetype="dashed", varSize=0.5, 
    drawVariance = TRUE, drawPsize = 0, 
    drawPlabel = FALSE, flagLimits=TRUE, 
    dirName="Direction", lineName="Experimental", ...)
{
  dots = list(...)
  p = list()
  gg = vario$getGgVec(idir,ivar,jvar) 
  hh = vario$getHhVec(idir,ivar,jvar)
  sw = vario$getSwVec(idir,ivar,jvar)
  df = data.frame(gg = gg, hh = hh, sw = sw)
  
  # Representing the Experimental variogram
  p = append(p, geom_line(data = df, 
  			mapping=aes(x=hh, y=gg, color=dirName, linetype=lineName), 
  			na.rm=TRUE, ...))
  
  # Representing the number of pairs (by size)
  if (drawPsize > 0)
  {
    p = append(p, geom_point(data = df, mapping=aes(x=hh, y=gg, size=sw), 
    		na.rm=TRUE, ...))
    p = append(p, list(labs(size = "Nb. pairs")))
  }
  else if (drawPsize < 0)
  {
    p = append(p, geom_point(data = df, mapping=aes(x=hh, y=gg), 
            na.rm=TRUE, ...))
  }
  
  # Representing the number of pairs (by label)
  if (drawPlabel)
    p = append(p, geom_text(data = df, mapping=aes(x=hh, y=gg, 
    		   label=as.character(sw)), na.rm=TRUE, ...))
  
  # Adding the vertical axis at X=0
  p = append(p, geom_vline(xintercept = 0., color='black', size=0.5))
  
  # Adding the horizontal axis at Y=0            
  p = append(p, geom_hline(yintercept = 0., color="black", size=0.5))
  
  # Drawing the variance-covariance reference line (optional)
  if (drawVariance)
    p = append(p, geom_hline(yintercept=vario$getVar(ivar,jvar), 
            color=varColor, linetype=varLinetype, size=varSize))
  
  # Tuning the bounds of graphics. This is optional in order to avoid multiple limit definitions
  if (flagLimits)
  {
    if (vario$drawOnlyPositiveX(ivar, jvar))
      p = append(p, plot.geometry(xlim = c(0, NA)))
    if (vario$drawOnlyPositiveY(ivar, jvar))
      p = append(p, plot.geometry(ylim = c(0, NA)))
  }
  p
}

#' Represent an elementary Model
#' @noRd
.modelElementary <- function(model, ivar=0, jvar=0, codir=NA,
    nh = 100, hmax = NA, asCov=FALSE, flagEnvelop = TRUE, 
    envColor='black', envLinetype="dashed", envSize=0.5, 
    dirName = "Direction1", lineName="Model",
    ...)
{
  dots = list(...)
  p = list()
  ndim = model$getDimensionNumber()
  
  # if hmax not specified = 3*maximum range of the model's basic structures
  if (.isNotDef(hmax))
  {
    hmax = 0
    for (icova in 1:model$getCovaNumber())
    {
      range_max = max(model$getCova(icova-1)$getRanges())
      if (3*range_max > hmax)
        hmax = 3*range_max
    }
  }
  # if hmax is still not calculated, take it empirically equal to 1.
  if (.isNotDef(hmax) || hmax == 0) hmax = 1.
  
  if (.isNotDef(codir))
  {
    codir = rep(0, ndim)
    codir[1] = 1
  }
  
  istart = 0
  for (icova in 1:model$getCovaNumber())
  {
    if (model$getCovName(icova-1) == 'Nugget Effect')
      istart = 1 # do not plot the first lag (h=0) for nugget effect (discontinuity)
  }
  
  # Calculating distances 
  eps = hmax / 1000.
  hh = seq(eps, hmax, length.out=nh)
  
  # Represent the Model
  mode = CovCalcMode()
  mode$setAsVario(! asCov)
  gg = model$sample(hh, ivar=ivar, jvar=jvar, codir=codir, mode=mode)
  df = data.frame(gg = gg[istart:nh], hh = hh[istart:nh])
  p = append(p, geom_line(data = df, 
  		     mapping=aes(x=hh, y=gg, color=dirName, linetype=lineName), 
             na.rm=TRUE, ...))
  
  # Represent the coregionalization envelop
  if (ivar != jvar && flagEnvelop)
  {
    gg = model$envelop(hh, ivar=ivar, jvar=jvar, isign=-1, codir=codir, mode=mode)
    df = data.frame(gg = gg[istart:nh], hh = hh[istart:nh])
    p = append(p, geom_line(data = df, mapping=aes(x=hh, y=gg), na.rm=TRUE, 
            color = envColor, linetype = envLinetype, size=envSize))
    
    gg = model$envelop(hh, ivar=ivar, jvar=jvar, isign=+1, codir=codir, mode=mode)
    df = data.frame(gg = gg[istart:nh], hh = hh[istart:nh])
    p = append(p, geom_line(data = df, mapping=aes(x=hh, y=gg), na.rm=TRUE, 
            color = envColor, linetype = envLinetype, size=envSize))
  }
  p
}

#' Represent an experimental variogram and overlay the Model calculated in same conditions
#' on a single figure
#' 
#' @param vario An object of the class Vario of gstlearn (optional)
#' @param model An object of the class Model of gstlearn (optional)
#' @param ivar Rank of the variable to be represented (-1 for all variables)
#' @param jvar Rank of the second variable to be represented in multivariate case (-1 for all variables
#' @param idir Rank of the direction to be represented (-1 for all directions)
#' @param nh Number of distance lags (used if 'vario' is not defined)
#' @param hmax Maximum distance (used if 'vario' is not defined)
#' @param drawPsize Represent variogram lags with a symbol proportional to the number of pairs
#' @param drawPlabel Represent variogram lags with the number of pairs displayed
#' @param asCov Represent the variogram as a covariance
#' @param drawVariance Represent statistical variance (or covariance)
#' @param flagEnvelop Represent the coregionalization envelop (multivariate case)
#' @param varioLinetype Linetype for representing the experimental variogram
#' @param modelLinetype Linetype for representing the Model
#' @param varColor Color for representing the variance (covariance)
#' @param varLinetype Linetype used for representing variance (covariance)
#' @param varSize Dimension used for representing variance (covariance)
#' @param envColor Color used for representing coregionalization envelop
#' @param envLinetype Linetype used for representing coregionalization envelop
#' @param envSize Size used for representing coregionalization envelop
#' @param drawVario Flag for representing the experimental variogram (used when 'vario' is defined)
#' @param flagLegend Flag for displaying the legend
#' @param ... Arguments passed to varioElementary() and modelElementary()
#' @return The ggplot object
plot.varmod <- function(vario=NA, model=NA, ivar=0, jvar=0, idir=-1,
    nh = 100, hmax = NA, drawPsize=-1, drawPlabel=FALSE, 
    asCov=FALSE, drawVariance = TRUE, flagEnvelop=TRUE, 
    varioLinetype = "dashed", modelLinetype = "solid",
    varColor='black', varLinetype="dashed", varSize=0.5, 
    envColor='black', envLinetype="dashed", envSize=0.5,
    drawVario=TRUE, flagLegend=FALSE, ...)
{
  dots = list(...)
  has_color = "color" %in% names(dots)
  has_codir = "codir" %in% names(dots)
  linetypes = c(varioLinetype, modelLinetype)
  
  p = list()
  ndir = 1
  if (! .isNotDef(vario)) ndir = vario$getDirectionNumber()
  nvar = 1
  if (! .isNotDef(vario)) nvar = vario$getVariableNumber()
  if (! .isNotDef(model)) nvar = model$getVariableNumber()
  cols = .getColors()
  
  idirUtil = .selectItemsInList(ndir, idir)
  ivarUtil = .selectItemsInList(nvar, ivar)
  jvarUtil = .selectItemsInList(nvar, jvar)
  ivarN = length(ivarUtil)
  jvarN = length(jvarUtil)
  
  if (.isNotDef(hmax))
  {
    if (! .isNotDef(vario))
      hmax = vario$getHmax(ivar, jvar, idir)
    else
      hhmax = 1
  }
  
  # Loop on the variables
  flag_allow_negative_X = FALSE
  flag_allow_negative_Y = FALSE
  
  # Allow redefining color and linetypes
  p <- append(p, list(new_scale("color")))
  p <- append(p, list(new_scale("linetype")))
 
  for (ivar in ivarUtil)
  {
    for (jvar in jvarUtil)
    {
      
      # Define the current plot
      for (idir in idirUtil)
      {
      	if (! .isNotDef(vario))
	      dirName = paste("Vario dir =", paste(round(vario$getCodirs(idir),3), collapse=' '))
 	    else
 	      dirName = paste("Direction :",idir)
    	
        # Plotting the experimental variogram
        if (! .isNotDef(vario) && drawVario)
        {
          dotloc = dots
          p = append(p, do.call(.varioElementary, c(list(vario=vario, 
                          ivar=ivar, jvar=jvar, idir=idir, 
                          varColor=varColor, varLinetype=varLinetype, varSize=varSize,
                          drawVariance=drawVariance, drawPsize=drawPsize, 
                          drawPlabel=drawPlabel, 
                          flagLimits=FALSE, dirName=dirName), dotloc)))
        }
        
        # Plotting the Model
        if (! .isNotDef(model))
        {
          dotloc = dots
          if (! .isNotDef(vario) && ! has_codir) dotloc$codir = vario$getCodirs(idir) 
          p = append(p, do.call(.modelElementary, c(list(model, ivar, jvar,  
                          nh = nh, hmax = hmax, asCov=asCov, 
                          flagEnvelop=flagEnvelop,
                          envColor = envColor, envLinetype = envLinetype, 
                          envSize=envSize, dirName=dirName), dotloc)))
        }
      } # End of loop on idir
             
      # Informing bound criterion
      if (! .isNotDef(vario))
      {
        if (! vario$drawOnlyPositiveX(ivar, jvar))
          flag_allow_negative_X = TRUE
        if (! vario$drawOnlyPositiveY(ivar, jvar))
          flag_allow_negative_Y = TRUE
      }
    } # End of loop on jvar
  } # End of loop on ivar
  
  # Adding some decoration
  p = append(p, plot.decoration(xlab = "Distance", ylab = "Variogram"))
  
  # Constructing the Legend

  p <- append(p, scale_linetype_manual(name="Types", values = linetypes))
  p <- append(p, scale_color_manual(name="Directions", values = cols))
  if (! flagLegend)
	p <- append(p, list(theme(legend.position='none')))
  
  # Tuning the global bounds of graphics
  lower_bound = NA
  if (! flag_allow_negative_X) lower_bound = 0
  p = append(p, plot.geometry(xlim = c(lower_bound, NA)))
  lower_bound = NA
  if (! flag_allow_negative_Y) lower_bound = 0
  p = append(p, plot.geometry(ylim = c(lower_bound, NA)))
  
  p
}

#' Arrange a set of figures for the multivariate representation
#' 
#' @param vario An object of the class Vario of gstlearn (optional)
#' @param model An object of the class Model of gstlearn (optional)
#' @param ivar Rank of the variable to be represented (-1 for all variables)
#' @param jvar Rank of the second variable to be represented in multivariate case (-1 for all variables
#' @param idir Rank of the direction to be represented (-1 for all directions)
#' @param nh Number of distance lags (used if 'vario' is not defined)
#' @param hmax Maximum distance (used if 'vario' is not defined)
#' @param drawPsize Represent variogram lags with a symbol proportional to the number of pairs
#' @param drawPlabel Represent variogram lags with the number of pairs displayed
#' @param asCov Represent the variogram as a covariance
#' @param drawVariance Represent statistical variance (or covariance)
#' @param flagEnvelop Represent the coregionalization envelop (multivariate case)
#' @param varioLinetype Linetype for representing the experimental variogram
#' @param modelLinetype Linetype for representing the Model
#' @param varColor Color for representing the variance (covariance)
#' @param varLinetype Linetype used for representing variance (covariance)
#' @param varSize Dimension used for representing variance (covariance)
#' @param envColor Color used for representing coregionalization envelop
#' @param envLinetype Linetype used for representing coregionalization envelop
#' @param envSize Size used for representing coregionalization envelop
#' @param drawVario Flag for representing the experimental variogram (used when 'vario' is defined)
#' @param flagLegend Flag for displaying the legend
#' @param ... Arguments passed to varioElementary() and modelElementary()
#' @return The ggplot object
multi.varmod <- function(vario, model=NA, ivar=-1, jvar=-1, idir=-1,
	nh = 100, hmax = NA, drawPsize=-1, drawPlabel=FALSE, 
    asCov=FALSE, drawVariance = TRUE, flagEnvelop=TRUE, 
    varioLinetype = "dashed", modelLinetype = "solid",
    varColor='black', varLinetype="dashed", varSize=0.5, 
    envColor='black', envLinetype="dashed", envSize=0.5,
    ...)
{
  nvar = vario$getVariableNumber()
  
  ivarUtil = .selectItemsInList(nvar, ivar)
  jvarUtil = .selectItemsInList(nvar, jvar)
  ivarN = length(ivarUtil)
  jvarN = length(jvarUtil)
  
  if (.isNotDef(hmax))
    hmax = vario$getHmax(ivar, jvar, idir)
  
  # Loop on the variables
  
  index = 0
  plot_lst <- vector("list", length = ivarN * jvarN)
  
  for (ivar in ivarUtil)
    for (jvar in jvarUtil)
    {
      
      # Define the current plot
      index = index + 1
      
      g = ggDefault()
      
      if (ivar < jvar)
      {
        g = g + theme(
            panel.background = element_rect(fill='transparent'), #transparent panel bg
            plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank(), #remove minor gridlines
            legend.background = element_rect(fill='transparent'), #transparent legend bg
            legend.box.background = element_rect(fill='transparent')) #transparent legend panel
      }
      else
      {
        g = g + plot.varmod(vario=vario, model=model, ivar=ivar, jvar=jvar, idir=-1,
            nh = nh, hmax = hmax, drawPsize=drawPsize, drawPlabel=drawPlabel, 
            asCov=asCov, drawVariance = drawVariance, flagEnvelop=flagEnvelop, 
            varioLinetype=varioLinetype, modelLinetype=modelLinetype,
            varColor=varColor, varLinetype=varLinetype, varSize=varSize, 
            envColor=envColor, envLinetype=envLinetype, envSize=envSize,
            ...)
      }
      plot_lst[[index]] <- g
    } 
  p = ggarrange(plotlist=plot_lst, nrow=ivarN, ncol = jvarN)
  p
}

#' Read a set of sample coordinates
#' @param db Db item from the gstlearn library
#' @param useSel Use of an optional selection (masking off samples)
#' @param posX Rank of the coordinate which will serve as first coordinate
#' @param posY Rank of the coordinate which will serve as the second coordinate
#' @return a Dataframe containing the 2-D coordinates
#' @noRd
.readPointCoor <- function(db, useSel=TRUE, posX=0, posY=1)
{
  if (db$getNDim() > 0) 
	x = db$getCoordinates(posX,useSel)
  if (db$getNDim() > 1)
    y = db$getCoordinates(posY,useSel)
  df = data.frame(x,y)
  df
}

#' Read the coordinates and one variable along a section of a Grid
#' @param dbgrid Grid data base from the gstlearn library
#' @param name Name of the target variable
#' @param useSel Use of an optional selection
#' @param posX rank of the coordinate which will serve as the first coordinate
#' @param posY rank of the coordinate which will serve as the second coordinate
#' @param corner A vector (same space dimension as 'dbgrid') which defines a pixel belonging to the extracted section
#' @return A dataframe containing the 2-D coordinates and the target variable
#' @note: setting useSel to FALSE enables having information for the whole grid (which remains a regular grid) even if a selection is defined.
#' @noRd
.readGridCoor <- function(dbgrid, name, useSel= FALSE, posX=0, posY=1, corner=NA)
{
  if (.isNotDef(corner))
	corner = rep(0, dbgrid$getNDim())
  
  if (dbgrid$getNDim() == 1)
  {
  	data = dbgrid$getColumn(name, useSel, FALSE)
  	x = dbgrid$getColumnByLocator(ELoc_X(), posX, FALSE, FALSE)
  	y = dbgrid$getColumnByLocator(ELoc_X(), posY, FALSE, FALSE)
  }
  else
  {
  	data = dbgrid$getOneSlice(name, posX, posY, corner, useSel)
  	nameX = dbgrid$getNameByLocator(ELoc_X(), posX)
  	x = dbgrid$getOneSlice(nameX, posX, posY, corner, FALSE)
  	nameY = dbgrid$getNameByLocator(ELoc_X(), posY)
  	y = dbgrid$getOneSlice(nameY, posX, posY, corner, FALSE)
  }
  	
  if (length(data) != length(x))
  {
    cat("Variable",name,"does not exist or does not have correction dimension\n")
    stop()
  }
  df = data.frame(x,y,data)
  df
}

#' Plotting a point data base where samples are displayed with different color and/or size
#' @param db Data Base containing the information to be displayed
#' @param nameColor Name of the variable to be represented in color
#' @param nameSize Name of the variable to be represented in proportional symbols
#' @param flagAbsSize Using the absolute value of the variable for graphic representation
#' @param flagCst Represent the location of the active samples only
#' @param useSel Use of the optional selection
#' @param asFactor Transform color variable into factor to use discrete palette
#' @param posX Rank of the coordinate used as the first coordinate
#' @param posY Rank of the coordinate used as the second coordinate
#' @param ... List of arguments passed to geom_point()
#' @return The description of the contents of the graphic layer
pointSymbol <- function(db, nameColor=NULL, nameSize=NULL,
    flagAbsSize = FALSE, flagCst=FALSE, useSel=TRUE, asFactor=FALSE, 
    posX=0, posY=1,
    ...) 
{ 
  # Creating the necessary data frame
  df = .readPointCoor(db, useSel, posX, posY)
  
  # Color of symbol
  colval = "constant"
  if (! is.null(nameColor)) {
    colval  = db$getColumn(nameColor, TRUE)
    if (asFactor) colval = factor(colval)
  }
  df["colval"] = colval 
  
  # Size of symbol
  sizval = "constant"
  if (! is.null(nameSize) && ! flagCst) {
    sizval  = db$getColumn(nameSize, TRUE)
    if (flagAbsSize) sizval = abs(sizval)
  }
  df["sizval"] = sizval
  
  if (flagCst)
 	 layer <- geom_point(data = df, mapping = aes(x=x, y=y), 
    		  na.rm=TRUE, ...)
  else
  	 layer <- geom_point(data = df, mapping = aes(x=x, y=y, color=colval, size=sizval), 
    		  na.rm=TRUE, ...)
  
  layer
}

#' Plotting a point data base where samples are displayed with a label
#' @param db Data Base containing the information to be displayed
#' @param name Name of the variable to be represented
#' @param digit Number of decimal digits
#' @param useSel Use of the optional selection
#' @param posX Rank of the coordinate used as the first coordinate
#' @param posY Rank of the coordinate used as the second coordinate
#' @param ... List of arguments passed to geom_text()
#' @return The description of the contents of the graphic layer
pointLabel <- function(db, name, digit=2, useSel=TRUE, posX=0, posY=1, ...) 
{  
  # Creating the necessary data frame
  df = .readPointCoor(db, useSel, posX, posY)
  
  # Label of symbols
  labval  = round(db$getColumn(name,TRUE),digit)
  df["labval"] = as.character(labval)
  
  layer <- geom_text_repel(data = df, mapping=aes(x=x, y=y, label=labval, 
           color=name), na.rm=TRUE, ...)
  
  layer
}

#' Plotting a point data base
#' @param db Data Base containing the information to be displayed
#' @param nameColor Name of the variable to be represented in color
#' @param nameSize Name of the variable to be represented in proportional symbols
#' @param nameLabel Name of the variable to be represented in literal manner
#' @param sizmin Minimum symbol size for proportional representation
#' @param sizmax Maximum symbol size for proportional representation
#' @param flagAbsSize Using the absolute value of the variable for graphic representation
#' @param flagCst Represent the location of the active samples only
#' @param palette Name of the reference color map
#' @param asFactor Transform the color variable into factor in order to use discrete palette
#' @param flagLegendColor Display the legend for Color representation
#' @param flagLegendSize Display the legend for Size representation
#' @param flagLegendLabel Display the legend for literal representation
#' @param legendNameColor Name of the Legend for color representation (set to 'nameColor' if not defined)
#' @param legendNameSize Name of the Legend for proportional representation (set to 'nameSize' if not defined)
#' @param legendNameLabel Name of the Legend for Literal representation (set to 'nameLabel' if not defined)
#' @param ... List of arguments passed to pointSymbol( ) and pointLabel() 
#' @return The ggplot object
plot.point <- function(db, nameColor=NULL, nameSize=NULL, nameLabel=NULL,
    sizmin=1, sizmax=5, flagAbsSize = FALSE, flagCst=FALSE, palette=NULL, asFactor=FALSE, 
    flagLegendColor=FALSE, flagLegendSize=FALSE, flagLegendLabel=FALSE, 
    legendNameColor=NULL, legendNameSize=NULL, legendNameLabel=NULL, 
    textColor="black", ...)
{ 
  p = list()
  title = ""
  
# If no variable is defined, use the default variable for Symbol(size) representation
# The default variable is the first Z-locator one, or the last variable in the file
  flagTitleDefault = FALSE
  if (is.null(nameColor) && is.null(nameSize) && is.null(nameLabel))
  {
  	nameSize = .getDefaultVariable(db)
    if (db$getLocNumber(ELoc_Z()) > 0)
      nameSize = db$getNameByLocator(ELoc_Z(),0)
    else 
    {
      # if no Z locator, choose the last field
      nameSize = db$getLastName()
      flagCst = TRUE
      flagTitleDefault = TRUE
    }
  }
  
  # Allow redefining color and linetypes
  p <- append(p, list(new_scale_color()))
  
  if (! is.null(nameColor) || ! is.null(nameSize))
  {
    p <- append(p, pointSymbol(db, nameColor=nameColor, nameSize=nameSize,
        flagAbsSize = flagAbsSize, flagCst=flagCst, asFactor=asFactor,
        ...))
    
    if (! is.null(nameSize) && ! flagCst)
    {
      p <- append(p, scale_size_continuous(range = c(sizmin, sizmax)))
    }
      
	# Palette definition (if defined)
	if (! is.null(palette))
	  p <- append(p, .scaleColorFill(palette, ...))
    
    # Set the default title
    if (! is.null(nameColor))
      title = paste(title, nameColor, "(color)", sep=" ")
    if (! is.null(nameSize))
      title = paste(title, nameSize, "(size)", sep=" ")
    
    # Set the Legend
    if (flagLegendColor)
    {
    	if (is.null(legendNameColor)) legendNameColor = nameColor
	    p <- append(p, list(labs(color = legendNameColor)))
	}
	else
	{
		p <- append(p, list(guides(color = "none")))
	}
	if (flagLegendSize && ! flagCst)
	{
		if (is.null(legendNameSize)) legendNameSize = nameSize
	    p <- append(p, list(labs(size = legendNameSize)))
	}
	else
	{
		p <- append(p, list(guides(size = "none")))
	}
  }
  
  if (! is.null(nameLabel))
  {
    p <- append(p, pointLabel(db, name=nameLabel, ...))
    p <- append(p, scale_color_manual(values = textColor))
    
    # Set the title              
    title = paste(title, nameLabel, sep=" ")
    
    # Set the legend
    if (flagLegendLabel)
    {
    	if (is.null(legendNameLabel)) legendNameLabel = nameLabel
	    p <- append(p, list(labs(color = legendNameLabel)))
	}
		else
	{
		p <- append(p, list(guides(color = "none")))
	}
  }

  
  # Decoration
  if (flagTitleDefault) title = "Sample Location"
  p <- append(p, plot.decoration(title = title))
    
  p
}

#' Represent the contents of a variable defined on a grid as an Image
#' @param dbgrid Grid data base from gstlearn
#' @param name Name of the variable to be represented
#' @param useSel Use of an optional selection
#' @param posX rank of the coordinate which will serve as the first coordinate
#' @param posY rank of the coordinate which will serve as the second coordinate
#' @param corner A vector (same space dimension as 'dbgrid') which defines a pixel belonging to the extracted section
#' @param ... Arguments passed to geom_tile() or geom_polygon()
#' @return The description of the contents of the figure
gridRaster <- function(dbgrid, name, useSel = TRUE, posX=0, posY=1, corner=NA, ...)
{
  # Reading the Grid information
  df = .readGridCoor(dbgrid, name, useSel, posX, posY, corner)
  
  # Define the contents
  if (dbgrid$getAngles()[1] == 0 && ! dbgrid$hasSingleBlock())
  {
    layer <- geom_tile(data = df, mapping=aes(x = x, y = y, fill = data), ...)
  }
  else
  {
    ids = seq(1, dbgrid$getNTotal())
    coords = dbgrid$getAllCellsEdges()
    positions = data.frame(id = rep(ids, each=4), x=coords[[1]], y=coords[[2]])
    values = data.frame(id = ids, value = df$data)
    df <- merge(values, positions, by = c("id"))
    layer <- geom_polygon(data = df, mapping=aes(x = x, y = y, fill = value, group = id), ...)
  }
  layer
}

#' Represent the contents of a variable defined on a grid with isovalues 
#' @param dbgrid Grid data base from gstlearn
#' @param name Name of the variable to be represented
#' @param useSel Use of an optional selection
#' @param posX rank of the coordinate which will serve as the first coordinate
#' @param posY rank of the coordinate which will serve as the second coordinate
#' @param corner A vector (same space dimension as 'dbgrid') which defines a pixel belonging to the extracted section
#' @param ... Arguments passed to geom_contour()
#' @return The description of the contents of the figure
gridContour <- function(dbgrid, name, useSel = TRUE, posX=0, posY=1, corner=NA, ...)
{
  # Reading the Grid information
  df = .readGridCoor(dbgrid, name, useSel, posX, posY, corner)
  
  layer <- geom_contour(data = df, mapping=aes(x = x, y = y, z = data), ...)
  
  layer
}

#' Define the variable used by default (when not explictly defined)
#' @param db Data base from gstlearn
#' @note The defaulted variable is the first one attached to the locator-Z (if any);
#' @note otherwise it is the last defined variable within 'db'.
#' @return Name of the defaulted variable
#' @noRd
.getDefaultVariable <- function(db)
{
	if (db$getLocNumber(ELoc_Z()) > 0)
		name = db$getNameByLocator(ELoc_Z(),0)
 	else
   	   # if no Z locator, choose the last field
  		name = db$getLastName()
  	name
}

#' Plotting a grid data base
#' @param dbgrid Grid Data Base containing the information to be displayed
#' @param nameRaster Name of the variable to be represented as an image
#' @param nameContour Name of the variable to be represented in isovalues
#' @param useSel Use of an optional selection
#' @param palette Name of the reference color map
#' @param naColor Color assigned to undefined samples
#' @param limits Bounds applied to the variable to be represented
#' @param flagLegendRaster Display the legend for grid representation as an image
#' @param flagLegendContour Display the legend for grid representation as contour lines
#' @param legendNameRaster Name of the Legend for representation as an image (set to 'nameRaster' if not defined)
#' @param legendNameContour Name of the Legend for representation as contour lines (set to 'nameContour' if not defined)
#' @param ... Arguments passed to gridRaster( ), gridContour(), .scaleColorFill()
#' @return The ggplot object
plot.grid <- function(dbgrid, nameRaster=NULL, nameContour=NULL,
    useSel = TRUE, palette=NULL, naColor = "white", limits = NULL, 
    flagLegendRaster=FALSE, flagLegendContour=FALSE,
    legendNameRaster=NULL, legendNameContour=NULL,
    ...)
{
  if (! dbgrid$isGrid())
  {
    cat("This function is restricted to Grid Db and cannot be used here\n")
    stop()
  }
  
  p = list()
  title = ""
  
  # If no variable is defined, use the default variable for Raster representation
  # The default variable is the first Z-locator one, or the last variable in the file
  if (is.null(nameRaster) && is.null(nameContour))
  {
  	nameRaster = .getDefaultVariable(dbgrid)
  }
  
  # Allow redefining color and linetypes
  p <- append(p, list(new_scale_color()))
  
  # Raster representation
  if (! is.null(nameRaster))
  {
    p <- append(p, gridRaster(dbgrid, name=nameRaster, useSel=useSel, ...))
    
    # Set the title
    title = paste(title,nameRaster)
    
    # Set the Legend
    if (flagLegendRaster)
    {
    if (is.null(legendNameRaster)) legendNameRaster = nameRaster
      p <- append(p, list(guides(fill = guide_colorbar(title=legendNameRaster, reverse=FALSE))))
    }
    else
      p <- append(p, list(theme(legend.position='none')))
  }
  
  # Contour representation
  if (! is.null(nameContour))
  {
    p = append(p, gridContour(dbgrid, name=nameContour, useSel=useSel, ...))
    
    # Set the title                    
    title = paste(title, nameContour, sep=" ")
    
    # Set the Legend
    if (flagLegendContour)
    {
      if (is.null(legendNameContour)) legendNameContour = nameContour
      p <- append(p, list(labs(contour = legendNameContour)))
    }
    else
      p <- append(p, list(guides(contour = "none")))
  }  
  
  # Palette definition (if defined)
  if (! is.null(palette))
	p <- append(p, .scaleColorFill(palette, naColor=naColor, limits=limits, ...))
  
  # Decoration
  p <- append(p, plot.decoration(title = title))
  
  p
}

#' Representing a Polygon
#' @param poly A Polygon object from the gstlearn library
#' @param flagTitle Plot the decoration attached to the figure
#' @param ... List of arguments passed to geom_polygon( )
#' @return The ggplot object
plot.polygon <- function(poly, flagTitle=FALSE, ...)
{
  dots = list(...)
  has_color = "color" %in% names(dots)
  
  p = list()
  npol = poly$getPolyElemNumber()
  cols = .getColors()
  
  dotloc = dots
  for (ipol in 1:npol)
  {
    if (! has_color) dotloc$color=cols[ipol+1]
    df = data.frame(x = poly$getX(ipol-1), y = poly$getY(ipol-1))
    p <- append(p, do.call(geom_polygon, c(list(data = df, mapping=aes(x=x,y=y)),  dotloc)))
  }
  
  # Decoration
  if (flagTitle)
    p <- append(p, plot.decoration(title = paste("Number of Polygons = ",npol)))
  
  p
}

#' Representing the histogram calculated on a variable
#' @param db A data base object from the gstlearn library
#' @param name Name of the target variable
#' @param useSel Use of an optional selection
#' @param ... List of arguments passed to geom_histogram( )
#' @return The ggplot object
plot.hist <- function(db, name, useSel=TRUE, ...)
{
  p = list()
  val  = db$getColumn(name, useSel)
  df = data.frame(val)
  
  p <- append(p, geom_histogram(data=df, mapping=aes(x=val), na.rm=TRUE, ...))
  p <- append(p, plot.decoration(title=name, xlab="Value", ylab="Count"))
  p
}

#' Representing the histogram calculated on an input vector of values
#' @param val The vector of input values
#' @param ... List of arguments passed to geom_histogram( )
#' @return The ggplot object
plot.histArray <- function(val, ...)
{
  p = list()  
  df = data.frame(val)
  
  p <- append(p, geom_histogram(data = df, mapping=aes(x=val), na.rm=TRUE, ...))
  p
}

#' Representing the 1-D polyline derived from the contents of a vector
#' @param data The vector containing the values to be represented
#' @param ... List of arguments passed to geom_line( )
#' @return The ggplot object
plot.curve <- function(data, ...)
{
  p = list()  
  absc = seq(1,length(data))
  df = data.frame(absc,data)
  
  p <- append(p, geom_line(data = df, mapping=aes(x=absc,y=data), na.rm=TRUE, ...))
  p
}

#' Representing the set of points
#' @param x Vector of first coordinate of the points to be visualized
#' @param y Vector of second coordinate of the points to be visualized
#' @param flagLine Option for joining the points to be displayed
#' @param flagPoint Option for representing the individual points to be displayed
#' @param ... List of arguments passed to geom_path() or geom_line()
#' @return The ggplot object
plot.XY <-function(x, y, flagLine=TRUE, flagPoint=FALSE, ...)
{
  if (length(y) != length(x))
  {
    cat("Arrays 'x' and 'y' should have same dimensions\n")
    stop()
  }
  
  p <- list()
  
  df = data.frame(x, y)
  
  if (flagLine)
    p <- append(p, geom_path(data = df, mapping=aes(x=x,y=y), na.rm=TRUE, ...))
  
  if (flagPoint)
    p <- append(p, geom_point(data = df, mapping=aes(x=x,y=y), na.rm=TRUE, ...))
  p
}

#' Representing the scatter plot defined as a 2-D grid of occurences
#' @param x Vector of first coordinate of the points to be visualized
#' @param y Vector of second coordinate of the points to be visualized
#' @param ... List of arguments passed to geom_bin2d()
#' @return The ggplot object
plot.hist2d <- function(x, y, ...)
{
  if (length(y) != length(x))
  {
    cat("Arrays 'x' and 'y' should have same dimensions\n")
    stop()
  }
  
  p <- list()
  df = data.frame(x, y)
  p = append(p, geom_bin2d(data = df, mapping = aes(x=x, y=y), na.rm=TRUE, ...))
  p = append(p, list(theme_bw()))
  p
}

#' Representing an anamorphosis
#' @param anam An experimental anamorphosis object from gstlearn library
#' @param ndisc Number of discretization points
#' @param aymin Minimum value along Gaussian axis
#' @param aymax Maximum value along Gaussian axis
#' @param ... List of arguments passed to plot.XY()
#' @return The ggplot object
plot.anam <- function(anam, ndisc=100, aymin=-10, aymax=10, ...)
{
  p = list()
  res = anam$sample(ndisc, aymin, aymax)
  
  p = append(p, plot.XY(res$getY(), res$getZ(), ...))
  p = append(p, plot.geometry(xlim=res$getAylim(), ylim=res$getAzlim()))
  p = append(p, plot.decoration(xlab = "Gaussian", ylab = "Raw"))
  p
}

#' Representing the scatter plot 
#' @param db1 A (first) data base from gstlearn library
#' @param namex Name of the variable (within 'db1') which will be displayed along the horizontal axis
#' @param namey Name of the variable (within 'db2') which will be displayed along the vertical axis
#' @param db2 A second data base from gstlearn. If not defined, it coincides with 'db1'
#' @param useSel Use of an optional selection (masking off samples)
#' @param asPoint Represent samples pointwise if TRUE, otherwise as a grid painted with occurrences
#' @param flagDiag Represent the diagonal of the plot
#' @param diagColor Color of the diagonal
#' @param diagLinetype Line type of the diagonal
#' @param flagRegr Represent the linear regression of Y|X 
#' @param regrColor Color of the linear regression of Y|X
#' @param regrLinetype Line type of the linear regression of Y|X
#' @param flagBiss Represent the first bisector (Y=X) 
#' @param bissColor Color of the first bisector (Y=X)
#' @param bissLinetype Line type of the first bisector (Y=X)
#' @param flagSameAxes Define the same bounds for horizontal and vertical axes
#' @param flagLegendRaster Show the legend when representing grid of occurrences (asPoint = FALSE)
#' @param legendNameRaster Name of the legend when representing grid of occurrences (asPoint = FALSE)
#' @param ... List of arguments passed to plot.XY() or plot.hist2d()
#' @return The ggplot object
plot.correlation <- function(db1, namex, namey, db2=NULL, useSel=TRUE,
    asPoint=FALSE, 
    flagDiag=FALSE, diagColor = "red", diagLinetype = "solid", 
    flagRegr=FALSE, regrColor = "blue", regrLinetype = "solid", 
    flagBiss=FALSE, bissColor = "green", bissLinetype = "solid", 
    flagSameAxes=FALSE, flagLegendRaster = FALSE, legendNameRaster = NULL,
    ...)
{
  if (is.null(db2)) db2 = db1
  x = db1$getColumn(namex, useSel)
  y = db2$getColumn(namey, useSel)
  
  p = list()
  if (asPoint)
    p = append(p, plot.XY(x, y, flagLine=FALSE, flagPoint=TRUE, ...))
  else
    p = append(p, plot.hist2d(x, y, ...))
  
  xmin = min(x, na.rm=TRUE)
  ymin = min(y, na.rm=TRUE)
  xmax = max(x, na.rm=TRUE)
  ymax = max(y, na.rm=TRUE)
  
  if (flagSameAxes)
  {
    xmin = ymin = min(xmin, ymin) 
    xmax = ymax = max(xmax, ymax)
  }
  
  if (flagDiag)
  {
    p <- append(p, geom_segment(aes(x=xmin,y=ymin,xend=xmax,yend=ymax),
        linetype = diagLinetype, color = diagColor, na.rm=TRUE))
  }
  
  if (flagBiss)
  {
    bmin = min(xmin, ymin)
    bmax = max(xmax, ymax)
    p <- append(p, geom_segment(aes(x=bmin,y=bmin,xend=bmax,yend=bmax),
        linetype = bissLinetype, color = bissColor, na.rm=TRUE))
  }
  
  if (flagRegr)
  {
    regr = regression(db2, namey, namex, flagCst=TRUE)
    if (regr$getNvar() > 0)
    {
      a = regr$getCoeff(0)
      b = regr$getCoeff(1)
      ymin = a + b * xmin
      ymax = a + b * xmax
      p <- append(p, geom_segment(aes(x=xmin,y=ymin,xend=xmax,yend=ymax),
          linetype = regrLinetype, color = regrColor, na.rm=TRUE))
    }
  }
  
  p = append(p, plot.decoration(xlab=namex, ylab=namey))
  
  # Set the Legend
  if (flagLegendRaster)
  {
    if (is.null(legendNameRaster)) legendNameRaster = "Count"
    p <- append(p, list(guides(fill = guide_colorbar(title=legendNameRaster, reverse=FALSE))))
  }
  else
  {
    p <- append(p, list(theme(legend.position='none')))
  }
  
  p 
}

#' Display a lithotype rule
#' @param rule A Rule object from gstlearn library
#' @param proportions The vector of facies proportions. When defined it is used to dimension the facies rectangles
#' @param maxG Maximum gaussian value (in absolute value)
#' @param ... List of arguments passed to geom_rect()
#' @return The ggplot object
plot.rule <- function(rule, proportions=NULL, maxG = 3., ...)
{
  p = list()
  nrect = rule$getFaciesNumber()
  if (! is.null(proportions)) 
    rule$setProportions(proportions)
  else
    rule$setProportions()
  cols = .getColors()
  
  df = data.frame(xmin=rep(0,nrect),xmax=rep(0,nrect),
      ymin=rep(0,nrect),ymax=rep(0,nrect),
      colors=cols[1:nrect])
  for (ifac in 1:nrect)
  {
    rect = rule$getThresh(ifac) # This function is 1-based
    df$xmin[ifac] = max(rect[1], -maxG)
    df$xmax[ifac] = min(rect[2], +maxG)
    df$ymin[ifac] = max(rect[3], -maxG)
    df$ymax[ifac] = min(rect[4], +maxG)
  }
  
  p = append(p, geom_rect(data = df, mapping=aes(xmin = xmin, xmax = xmax, 
          ymin = ymin, ymax = ymax, fill = colors), na.rm=TRUE, ...))
  p = append(p, plot.geometry(xlim=c(-maxG,+maxG), ylim=c(-maxG,+maxG)))
  p
}

#' Represent a meshing
#' @param mesh A mesh object from gstlearn library
#' @param flagFace Assign a color to each polygon of the meshing
#' @param flagApex Draw the rank of each apex of the meshing
#' @param rankMeshMax Rank of the maximum mesh to be represented (if > 0)
#' @param ... List of arguments passed to geom_polygon(), geom_path() or geom_point()
#' @return The ggplot object
plot.mesh <- function(mesh, flagFace=FALSE, flagApex=FALSE, rankMeshMax= -1, ...)
{
  p = list()
  
  nmesh = mesh$getNMeshes()
  nmax = nmesh
  if (rankMeshMax > 0) nmax = rankMeshMax
  for (imesh in 1:nmax)
  {
    x = mesh$getCoordinatesPerMesh(imesh-1, 0, TRUE)
    y = mesh$getCoordinatesPerMesh(imesh-1, 1, TRUE)
    df = data.frame(x, y)
    if (flagFace)
      p <- append(p, geom_polygon(data = df, mapping=aes(x=x,y=y), ...))
    else
      p <- append(p, geom_path(data = df, mapping=aes(x=x,y=y), ...))
    if (flagApex)
      p <- append(p, geom_point(data = df, mapping=aes(x=x, y=y), ...))
  }  
  p
}

#' Represent the neighborhood
#' @param neigh A Neigh object from the gstlearn library
#' @param grid The target Data base from the gstlearn library
#' @param node Rank of the target 
#' @param flagCell Represent the target as the corresponding cell
#' @param flagZoom Zoom to the extension of the neighborhood
#' @param ... Arguments passed to plot.XY()
#' @return The ggplot object
plot.neigh <- function(neigh, grid, node=0, flagCell=FALSE, flagZoom=FALSE, ...)
{
	p = list()

    # Identify target location
    target = grid$getSampleCoordinates(node)
    
    # Represent the target location
    p = append(p, plot.XY(target[1], target[2], flagLine=FALSE, flagPoint=TRUE, ...))
    
    # Represent the edge of the target (if block)
    edges = grid$getCellEdges(node)
    p = append(p, plot.XY(edges[[1]], edges[[2]], ...))
    
    # Represent the Neighborhood Ellipsoid
    if (neigh$getType()$getValue() == ENeigh_MOVING()$getValue())
    {
    	edges = neigh$getEllipsoid(target)
    	p = append(p, plot.XY(edges[[1]], edges[[2]], ...))
    }
    
    # Represent the Angular sectors
    if (neigh$getFlagSector())
    {
    	segments = neigh$getSectors(target)
        nseg = length(segments[[1]])
        cx = numeric(2)
        cx[1] = target[1]
        cy = numeric(2)
        cy[1] = target[2]
	    for (iseg in 1:nseg)
	    {
	    	cx[2] = segments[[1]][iseg]
	    	cy[2] = segments[[2]][iseg]
       		p = append(p, plot.XY(cx, cy, flagLine=TRUE, flagPoint=FALSE, ...))
 		}
    }
        
    # Zoom to the Maximum radius circle (optional)
    if (flagZoom)
    {
		limits = neigh$getZoomLimits(target)
		p <- p + plot.geometry(xlim=limits[[1]], ylim=limits[[2]]) 
    }
    p
}

#' Represent the non-stationary parameters of a Model on a Grid
#' @param model A Model object (carrying non-stationary parameters) from gstlearn
#' @param dbgrid A Grid data base used for representation
#' @param useSel Use of an optional selection (masking off samples)
#' @param icov Rank of the covariance used for display
#' @param color Color used for graphic representation
#' @param flagOrtho Defines the long_axis of the anisotropy with respect to angle
#' @param scale Size given to the arraws
#' @param ... Arguments passed to geom_segment()
#' @return The ggplot object to geom_segment()
plot.modelOnGrid <- function(model, dbgrid, useSel=TRUE, icov=0, color='black', 
	flagOrtho=TRUE, scale=40, ...)
{
    # Extracting coordinates
    tabx = dbgrid$getCoordinates(0,useSel)
    taby = dbgrid$getCoordinates(1,useSel)
    
    # Process the non-stationarity
    db_model_nostat(dbgrid, model, icov)
    tabR1 = dbgrid$getColumn("Nostat.Range-1", useSel)
    tabR2 = dbgrid$getColumn("Nostat.Range-2", useSel)
    tabA  = dbgrid$getColumn("Nostat.Angle-1", useSel)
    if (flagOrtho) tabA = 90 + tabA
    tabA = tabA * pi / 180.
    
    tabdx = (tabR1 * cos(tabA) - tabR2 * sin(tabA)) * scale
    tabdy = (tabR1 * sin(tabA) + tabR2 * cos(tabA)) * scale
    data = data.frame(x = tabx, y = taby, dx=tabdx, dy=tabdy)
#    ax.quiver(tabx, taby, tabR2, tabR2, angles=tabA, color=color, **kwargs)
    
  	p = ggplot(data = data, aes(x = x, y = y)) + 
 	   geom_point(size = 1) + 
 	   geom_segment(aes(xend = x + dx, yend = y + dy),
                 arrow = arrow(length = unit(0.1, "cm")), ...)
    
	p
}

# The following code has been moved in rgstlearn.i (to prevent roxygen from crashing
#setMethod("plot", signature(x="_p_AMesh"), function(x,y=missing,...)   plot.mesh(x,...))
#setMethod("plot", signature(x="_p_DbGrid"), function(x,y="missing",...)  plot.grid(x,...))

#setMethod("plot", signature(x="_p_Db"), function(x,y=missing,...) plot.point(x,...))
#setMethod("plot", signature(x="_p_Polygons"), function(x,y=missing,...) plot.polygon(x,...))

#setMethod("plot", signature(x="_p_Vario"), function(x,y=missing,...) plot.vario(x,...))
#setMethod("plot", signature(x="_p_Model"), function(x,y="missing",...) plot.model(x,...))

#setMethod("plot", signature(x="_p_Rule"), function(x,y="missing",...) plot.rule(x,...))
#setMethod("plot", signature(x="_p_AAnam"), function(x,y="missing",...) plot.anam(x,...))
