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
# plots easily 
#

#' Define the global values in the given environment position (search)
#
#' @param pos Position in the list of packages
plot.initialize <- function(pos=1) 
{
  assign("plot.default_dims", list(c(8,8), c(8,8)), pos=pos)
  assign("plot.default_xlim", list(c(NA,NA), c(NA,NA)), pos=pos)
  assign("plot.default_ylim", list(c(NA,NA), c(NA,NA)), pos=pos)
  assign("plot.default_asp",  c(0, 1), pos=pos)
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
  plot.setDefaultInternal(2, dims=dims, xlim=xlim, ylim=ylim, asp=asp)
}

#' Set the default values for all subsequent non-Geographical figures
#'
#' @param dims Vector giving the dimensions of the figure
#' @param xlim Bounds of the figure along the horizontal axis (when left to NA, it will be adjusted to the figure contents)
#' @param ylim Bounds of the figure along the vertical axis (when left to NA, it will be adjusted to the figure contents)
#' @param asp Aspect ratio Y/X
plot.setDefault <- function(dims=NA, xlim=NA, ylim=NA, asp=NA)
{
  plot.setDefaultInternal(1, dims=dims, xlim=xlim, ylim=ylim, asp=asp)
}

#' Set the default values for all subsequent Geographical figures
#' @param mode 1 for Geographical and 2 for non-Geographical parameteres
#' @param dims Dimensions of the figures
#' @param xlim Bounds along the horizontal axis
#' @param ylim Bounds along the vertical axis
#' @param asp Aspect ratio Y/X
#' @noRd
plot.setDefaultInternal <- function(mode=1, dims=NA, xlim=NA, ylim=NA, asp=NA)
{
  if (!.isNotDef(dims))
    plot.default_dims[[mode]] = dims
  if (!.isNotDef(xlim))
    plot.default_xlim[[mode]] = xlim
  if (!.isNotDef(ylim))
    plot.default_ylim[[mode]] = ylim    
  if (!.isNotDef(asp))
    plot.default_asp[[mode]] = asp
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
    
    if (!.isNotDef(plot.default_dims[[mode]]))
      cat("- Figure dimensions =", plot.default_dims[[mode]],"\n")
    else
      cat("- Figure dimensions (not defined)\n")
    
    if (!.isNotDef(plot.default_xlim[[mode]]))
      cat("- Limits along X =",plot.default_xlim[[mode]],"\n")
    else
      cat("- Limits along X (not defined)\n")
    
    if (!.isNotDef(plot.default_ylim[[mode]]))
      cat("- Limits along Y =",plot.default_ylim[[mode]],"\n")
    else
      cat("- Limits along Y (not defined)\n")
    
    if (plot.default_asp[mode] != 0)
      cat("- Aspect =",plot.default_asp[mode],"\n")
    else
      cat("- Aspect (automatic)\n")
  }  
}

#' Define the color map
#' @param palette Reference palette used for defining the current color map
#' @noRd
.scaleColorFill <- function(palette, ...)
{
  rcb <- c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", "PuBu",
      "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd",
      "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3",
      "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")
  v <- c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo",
      "A", "B", "C", "D", "E", "F", "G", "H")
  rcb_num <- 1:18
  
  aes_list = c("colour", "fill")
  if(length(palette) == 1) {
    if (any(palette == rcb) | any(palette == rcb_num)) {
      layer = scale_color_distiller(palette= palette, aesthetics=aes_list, ...)
    } else if(any(palette == v)) {
     layer = scale_color_viridis_c(option= palette, aesthetics=aes_list, ...)
    } 
  } else if(length(palette) == 2) {
    low = palette[1]
    high = palette[2]
    layer = scale_color_gradient(low= low, high= high, aesthetics=aes_list, ...)
  } else if(length(palette) == 3) {
    low = palette[1]
    mid = palette[2]
    high = palette[3]
    layer = scale_color_gradient2(low= low, mid= mid, high= high, aesthetics=aes_list, ...)
  } else {
    layer = scale_colour_manual(values= palette, aesthetics=aes_list, ...)
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
#' @param flag_suppress_warnings TRUE to suppress informational warnings
ggPrint <- function(p, flag_suppress_warnings = TRUE)
{
  if (flag_suppress_warnings)
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
 	locdims = plot.default_dims[[mode]]
  else
  	locdims = figsize
   
  p <- p + plot.geometry(dims=locdims, 
                         xlim=plot.default_xlim[[mode]], 
                         ylim=plot.default_ylim[[mode]], 
                         asp=plot.default_asp[mode])
  p
}

#' Initiate a non-geographical display (using the default values for parameters)
#' @return The initiated ggplot object
#' @note Use printDefault() to visualize them and setDefault() to modify them
ggDefault <- function(figsize=NA)
{
  p <- ggplot()
  mode = 1
  
  if (.isNotDef(figsize))
 	locdims = plot.default_dims[[mode]]
  else
  	locdims = figsize
  
  p <- p + plot.geometry(dims=locdims, 
                         xlim=plot.default_xlim[[mode]], 
                         ylim=plot.default_ylim[[mode]], 
                         asp=plot.default_asp[mode])
  p
}

#' Check if the argument can be considered as an array (with possibly required dimensions)
#' @param arg Input argument
#' @param ndim Required dimension for the input argument (no check is performed if NA)
#' @noRd
.isArray <- function(arg, ndim=NA)
{
  if (length(arg) <= 1) return (FALSE)
  
  if (!.isNotDef(ndim) && length(arg) != ndim) return (FALSE)
  
  TRUE
}

#' Draw the decoration of a figure (title, axis labels, ...)
#'
#' @param xlab Label along the horizontal axis
#' @param ymab Label along the vertical axis
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
  
  if (length(xlim) > 1)
  {
    if (.isArray(xlim, 2))
    {
      p <- c(p, scale_x_continuous(limits=xlim, expand=expand))
    }
    else
      cat("'xlim' should be [a,b] or [NA,b] or [a,NA]. Ignored\n")
  }
  
  if (length(ylim) > 1)
  {
    if (.isArray(ylim, 2))
    {
      p <- c(p, scale_y_continuous(limits=ylim, expand=expand))
    }
    else
      cat("'ylim' should be [a,b] or [NA,b] or [a,NA]. Ignored\n")
  }
  
  if (!.isNotDef(asp))
  {     
    if (asp != 0)
      p = c(p, coord_fixed(asp))
    else
      p = c(p, list(theme(aspect.ratio = 1)))
  }
  p
}

#' Function for representing a Model
#' @param model An object of class Model from gstlearn
#' @param ivar Rank of the variable to be represented (-1 for all variables)
#' @param jvar Rank of the second variable to be represented in multivariate case (-1 for all variables
#' @param vario An object of class Vario of gstlearn (optional)
#' @param idir Rank of the direction
#'
#' @notes If 'vario' is defined, the calculation direction is given by the definition of direction 'idir' within 'vario' 
#' @return The ggplot object
plot.model <- function(model, ivar=0, jvar=0, vario=NA, idir=0, ...)
{
  p = list()
  p = c(p, plot.varmod(vario=vario, model=model, ivar=ivar, jvar=jvar, idir=idir,
          draw_vario=FALSE, ...))
  p = c(p, plot.decoration(xlab = "Distance", ylab = "Variogram"))
  p  
}

#' Function for representing the Experimental Variogram
#' @param vario An object of the class Vario of gstlearn
#' @param ivar Rank of the variable to be represented 
#' @param jvar Rank of the second variable to be represented in multivariate case 
#' @param idir Rank of the direction to be represented 
#' @return The ggplot object
plot.vario <- function(vario, ivar=0, jvar=0, idir=0,...)
{
  p = list()
  p = c(p, plot.varmod(vario=vario, ivar=ivar, jvar=jvar, idir=idir,...))
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
    linetype="solid",
    var_color='black', var_linetype="dashed", var_size=0.5, 
    draw_variance = TRUE, draw_post=TRUE, draw_psize = 0, 
    draw_plabel = FALSE, label=NULL, flagLimits=TRUE, ...)
{
  dots = list(...)
  p = list()
  gg = vario$getGgVec(idir,ivar,jvar) 
  hh = vario$getHhVec(idir,ivar,jvar)
  sw = vario$getSwVec(idir,ivar,jvar)
  df = data.frame(gg = gg, hh = hh, sw = sw)
  
# Representing the Experimental variogram
#    args <- formals(geom_line)
#    args[which(names(args) %in% names(dots))] <- dots[na.omit(match(names(args), names(dots)))]
#    print(args)
#    layer <- do.call(geom_line, c(list(data = df, mapping=aes(x=hh, y=gg), na.rm=TRUE),args = args))
  p = c(p, geom_line(data = df, mapping=aes(x=hh, y=gg), na.rm=TRUE, 
          linetype=linetype, ...))
  
  # Representing the number of pairs (by size)
  if (draw_psize > 0)
  {
    p = c(p, geom_point(data = df, mapping=aes(x=hh, y=gg, size=sw), na.rm=TRUE, ...))
    p = c(p, list(labs(size = "Nb. pairs")))
  }
  else if (draw_psize < 0)
  {
    p = c(p, geom_point(data = df, mapping=aes(x=hh, y=gg), na.rm=TRUE, ...))
  }
  
  # Representing the number of pairs (by label)
  if (draw_plabel)
    p = c(p, geom_text(data = df, mapping=aes(x=hh, y=gg, label=as.character(sw)), 
            na.rm=TRUE, ...))
  
  # Constructing the label for Legend
  if (length(label) <= 0)
    label = paste("vario dir=", paste(vario$getCodirs(idir), collapse=' '))
  
  # Adding the vertical axis at X=0
  p = c(p, geom_vline(xintercept = 0., color='black', size=0.5))
  
  # Adding the horizontal axis at Y=0            
  p = c(p, geom_hline(yintercept = 0., color="black", size=0.5))
  
  # Drawing the variance-covariance reference line (optional)
  if (draw_variance)
    p = c(p, geom_hline(yintercept=vario$getVar(ivar,jvar), 
            color=var_color, linetype=var_linetype, size=var_size))
  
  # Tuning the bounds of graphics. This is optional in order to avoid multiple limit definitions
  if (flagLimits)
  {
    if (vario$drawOnlyPositiveX(ivar, jvar))
      p = c(p, plot.geometry(xlim = c(0, NA)))
    if (vario$drawOnlyPositiveY(ivar, jvar))
      p = c(p, plot.geometry(ylim = c(0, NA)))
  }
  p
}

#' Represent an elementary Model
#' @noRd
.modelElementary <- function(model, ivar=0, jvar=0, codir=NA,
    nh = 100, hmax = NA, asCov=FALSE, flag_envelop = TRUE, 
    env_color='black', env_linetype="dashed", env_size=0.5,
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
  p = c(p, geom_line(data = df, mapping=aes(x=hh, y=gg), na.rm=TRUE, ...))
  
  # Represent the coregionalization envelop
  if (ivar != jvar && flag_envelop)
  {
    gg = model$envelop(hh, ivar=ivar, jvar=jvar, isign=-1, codir=codir, mode=mode)
    df = data.frame(gg = gg[istart:nh], hh = hh[istart:nh])
    p = c(p, geom_line(data = df, mapping=aes(x=hh, y=gg), na.rm=TRUE, 
            color = env_color, linetype = env_linetype, size=env_size))
    
    gg = model$envelop(hh, ivar=ivar, jvar=jvar, isign=+1, codir=codir, mode=mode)
    df = data.frame(gg = gg[istart:nh], hh = hh[istart:nh])
    p = c(p, geom_line(data = df, mapping=aes(x=hh, y=gg), na.rm=TRUE, 
            color = env_color, linetype = env_linetype, size=env_size))
  }
  p
}

#' Represent an experimental variogram and overlay the Model calculated in same conditions
#' 
#' @param vario An object of the class Vario of gstlearn (optional)
#' @param model An object of the class Model of gstlearn (optional)
#' @param ivar Rank of the variable to be represented (-1 for all variables)
#' @param jvar Rank of the second variable to be represented in multivariate case (-1 for all variables
#' @param idir Rank of the direction to be represented (-1 for all directions)
#' @param nh Number of distance lags (used if 'vario' is not defined)
#' @param hmax Maximum distance (used if 'vario' is not defined)
#' @param draw_psize Represent variogram lags with a symbol proportional to the number of pairs
#' @param draw_plabel Represent variogram lags with the number of pairs displayed
#' @param asCov Represent the variogram as a covariance
#' @param draw_variance Represent statistical variance (or covariance)
#' @param flag_envelop Represent the coregionalization envelop (multivariate case)
#' @param var_color Color for representing the variance (covariance)
#' @param var_linetype Linetype used for representing variance (covariance)
#' @param var_size Dimension used for representing variance (covariance)
#' @param env_color Color used for representing coregionalization envelop
#' @param env_linetype Linetype used for representing coregionalization envelop
#' @param env_size Size used for representing coregionalization envelop
#' 0param draw_vario Flag for representing the experimental variogram (used when 'vario' is defined)
#' @param label Label defined for representation of the Legend
#' @return The ggplot object
plot.varmod <- function(vario=NA, model=NA, ivar=-1, jvar=-1, idir=-1,
    nh = 100, hmax = NA, draw_psize=-1, draw_plabel=FALSE, 
    asCov=FALSE, draw_variance = TRUE, flag_envelop=TRUE, 
    var_color='black', var_linetype="dashed", var_size=0.5, 
    env_color='black', env_linetype="dashed", env_size=0.5,
    draw_vario=TRUE, label=NULL, ...)
{
  dots = list(...)
  has_color = "color" %in% names(dots)
  has_linetype = "linetype" %in% names(dots)
  has_codir = "codir" %in% names(dots)
  
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
  
  for (ivar in ivarUtil)
    for (jvar in jvarUtil)
    {
      
      # Define the current plot
      for (idir in idirUtil)
      {
        # Plotting the experimental variogram
        if (! .isNotDef(vario) && draw_vario)
        {
          dotloc = dots
          if (! has_color) dotloc$color=cols[idir+1]
          if (! has_linetype) dotloc$linetype = "dashed"
          p = c(p, do.call(.varioElementary, c(list(vario=vario, ivar=ivar, jvar=jvar, idir=idir, 
                          var_color=var_color, var_linetype=var_linetype, var_size=var_size,
                          draw_variance=draw_variance, draw_psize=draw_psize, 
                          draw_plabel=draw_plabel, label=label, flagLimits=FALSE), dotloc)))
        }
        
        # Plotting the Model
        if (! .isNotDef(model))
        {
          dotloc = dots
          if (! has_color) dotloc$color=cols[idir+1]
          if (! has_linetype) dotloc$linetype = "solid"
          if (! .isNotDef(vario) && ! has_codir) dotloc$codir = vario$getCodirs(idir) 
          p = c(p, do.call(.modelElementary, c(list(model, ivar, jvar,  
                          nh = nh, hmax = hmax, asCov=asCov, flag_envelop=flag_envelop,
                          env_color = env_color, env_linetype = env_linetype, 
                          env_size=env_size), dotloc)))
        }
        
        # Adding some decoration
        p = c(p, plot.decoration(xlab = "Distance", ylab = "Variogram"))
      }
      
      # Informing bound criterion
      if (! .isNotDef(vario))
      {
        if (! vario$drawOnlyPositiveX(ivar, jvar))
          flag_allow_negative_X = TRUE
        if (! vario$drawOnlyPositiveY(ivar, jvar))
          flag_allow_negative_Y = TRUE
      }
    } 
  
  # Tuning the global bounds of graphics
  lower_bound = NA
  if (! flag_allow_negative_X) lower_bound = 0
  p = c(p, plot.geometry(xlim = c(lower_bound, NA)))
  lower_bound = NA
  if (! flag_allow_negative_Y) lower_bound = 0
  p = c(p, plot.geometry(ylim = c(lower_bound, NA)))
  
  p
}

#' Arrange a set of figures for the multivariate representation
#' Same arguments as in plot.varmod()
#' @return The ggplot object
multi.varmod <- function(vario, model=NA, ivar=-1, jvar=-1, idir=-1,
	nh = 100, hmax = NA, draw_psize=-1, draw_plabel=FALSE, 
    asCov=FALSE, draw_variance = TRUE, flag_envelop=TRUE, 
    var_color='black', var_linetype="dashed", var_size=0.5, 
    env_color='black', env_linetype="dashed", env_size=0.5,
    label=NULL, ...)
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
            nh = nh, hmax = hmax, draw_psize=draw_psize, draw_plabel=draw_plabel, 
            asCov=asCov, draw_variance = draw_variance, flag_envelop=flag_envelop, 
            var_color=var_color, var_linetype=var_linetype, var_size=var_size, 
            env_color=env_color, env_linetype=env_linetype, env_size=env_size,
            label=label, ...)
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
#' @param name_color Name of the variable to be represented in color
#' @param name_size Name of the variable to be represented in proportional symbols
#' @param flagAbsSize Using the absolute value of the variable for graphic representation
#' @param flagCst Represent the location of the active samples only
#' @param useSel Use of the optional selection
#' @param asFactor Transform color variable into factor to use discrete palette
#' @param posX Rank of the coordinate used as the first coordinate
#' @param posY Rank of the coordinate used as the second coordinate
#' @param ... List of arguments passed to geom_point()
#' @return The description of the contents of the graphic layer
pointSymbol <- function(db, name_color=NULL, name_size=NULL,
    flagAbsSize = FALSE, flagCst=FALSE, useSel=TRUE, asFactor=FALSE, posX=0, posY=1,
    ...) 
{ 
  # Creating the necessary data frame
  df = .readPointCoor(db, useSel, posX, posY)
  
  # Color of symbol
  colval = NULL
  if (! is.null(name_color)) {
    colval  = db$getColumn(name_color, TRUE)
    if (asFactor) colval = factor(colval)
  }
  df["colval"] = colval 
  
  # Size of symbol
  sizval = NULL
  if (! is.null(name_size)) {
    if (! flagCst)
    {
      sizval  = db$getColumn(name_size, TRUE)
      if (flagAbsSize) sizval = abs(sizval)
    }
  }
  df["sizval"] = sizval
  
  layer <- geom_point(data = df, mapping = aes(x=x, y=y, color=colval, size=sizval), 
      na.rm=TRUE, ...)
  
  layer
}

#' Plotting a point data base where samples are displayed with a label
#' @param db Data Base containing the information to be displayed
#' @param name Name of the variable to be represented
#' @param digits Number of decimal digits
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
  
  layer <- geom_text(data = df, mapping=aes(x=x, y=y, label=labval), 
      na.rm=TRUE, ...)
  
  layer
}

#' Plotting a point data base
#' @param db Data Base containing the information to be displayed
#' @param name_color Name of the variable to be represented in color
#' @param name_size Name of the variable to be represented in proportional symbols
#' @param name_label Name of the variable to be represented in literal manner
#' @param sizmin Minimum symbol size for proportional representation
#' @param sizmax Maximum symbol size for proportional representation
#' @param flagAbsSize Using the absolute value of the variable for graphic representation
#' @param flagCst Represent the location of the active samples only
#' @param palette Name of the reference color map
#' @param asFactor Transform the color variable into factor in order to use discrete palette
#' @param show.legend.color Display the legend for Color representation
#' @param show.legend.size Display the legend for Size representation
#' @param show.legend.label Display the legend for literal representation
#' @param legend.name.color Name of the Legend for color representation (set to 'name_color' if not defined)
#' @param legend.name.size Name of the Legend for proportional representation (set to 'name_size' if not defined)
#' @param legend.name.label Name of the Legend for Literal representation (set to 'name_literam' if not defined)
#' @param ... List of arguments passed to pointSymbol( ) and pointLabel() 
#' @return The ggplot object
plot.point <- function(db, name_color=NULL, name_size=NULL, name_label=NULL,
    sizmin=1, sizmax=5, flagAbsSize = FALSE, flagCst=FALSE, palette=NULL, asFactor=FALSE, 
    show.legend.color=FALSE, show.legend.size=FALSE, show.legend.label=FALSE, 
    legend.name.color=NULL, legend.name.size=NULL, legend.name.label=NULL, ...)
{ 
  p = list()
  title = ""
  
# If no variable is defined, use the default variable for Symbol(size) representation
# The default variable is the first Z-locator one, or the last variable in the file
  flagTitleDefault = FALSE
  if (is.null(name_color) && is.null(name_size) && is.null(name_label))
  {
  	name_size = .getDefaultVariable(db)
    if (db$getLocNumber(ELoc_Z()) > 0)
      name_size = db$getNameByLocator(ELoc_Z(),0)
    else 
    {
      # if no Z locator, choose the last field
      name_size = db$getLastName()
      flagCst = TRUE
      flagTitleDefault = TRUE
    }
  }
  
  if (! is.null(name_color) || ! is.null(name_size))
  {
    p <- c(p, pointSymbol(db, name_color=name_color, name_size=name_size,
        flagAbsSize = flagAbsSize, flagCst=flagCst, asFactor=asFactor,
        ...))
    
    if (! is.null(name_size) && ! flagCst)
      p <- c(p, scale_size_continuous(range = c(sizmin, sizmax)))
      
	# Palette definition
	if (! is.null(palette))
    	p <- c(p, .scaleColorFill(palette, ...))
    
    # Set the default title
    if (! is.null(name_color))
      title = paste(title, name_color)
    if (! is.null(name_size))
      title = paste(title, name_size, sep=" ")
    
    # Set the Legend
    if (show.legend.color)
    {
    	if (is.null(legend.name.color)) legend.name.color = name_color
	    p <- c(p, list(labs(color = legend.name.color)))
	}
	if (show.legend.size)
	{
		if (is.null(legend.name.size)) legend.name.size = name_size
	    p <- c(p, list(labs(size = legend.name.size)))
	}
  }
  
  if (! is.null(name_label))
  {
    p <- c(p, pointLabel(db, ...))
    
    # Set the title              
    title = paste(title, name_label, sep=" ")
    
    # Set the legend
    if (show.legend.label)
    {
    	if (is.null(legend.name.label)) legend.name.label = name_label
	    p <- c(p, list(labs(label = legend.name.label)))
	}
  }
  
  # Decoration
  if (flagTitleDefault) title = "Sample Location"
  p <- c(p, plot.decoration(title = title))
  
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
#' @param name_raster Name of the variable to be represented as an image
#' @param name_contour Name of the variable to be represented in isovalues
#' @param useSel Use of an optional selection
#' @param palette Name of the reference color map
#' @param na.value Color assigned to undefined samples
#' @param limits Bounds applied to the variable to be represented
#' @param show.legend.raster Display the legend for grid representation as an image
#' @param show.legend.contour Display the legend for grid representation as contour lines
#' @param legend.name.raster Name of the Legend for representation as an image (set to 'name_raster' if not defined)
#' @param legend.name.contour Name of the Legend for representation as contour lines (set to 'name_contour' if not defined)
#' @param ... List of arguments passed to gridRaster( ), gridContour() and .scaleColorFill()
#' @return The ggplot object
plot.grid <- function(dbgrid, name_raster=NULL, name_contour=NULL,
    useSel = TRUE, palette=NULL, na.value = "white", limits = NULL, 
    show.legend.raster=FALSE, show.legend.contour=FALSE,
    legend.name.raster=NULL, legend.name.contour=NULL,
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
  if (is.null(name_raster) && is.null(name_contour))
  {
  	name_raster = .getDefaultVariable(dbgrid)
  }
  
  # Raster representation
  if (! is.null(name_raster))
  {
    p <- c(p, gridRaster(dbgrid, name=name_raster, useSel=useSel, ...))
    
    # Set the title
    title = paste(title,name_raster)
    
    # Set the Legend
    if (show.legend.raster)
    {
    if (is.null(legend.name.raster)) legend.name.raster = name_raster
      p <- c(p, list(guides(fill = guide_colorbar(title=legend.name.raster, reverse=FALSE))))
    }
    else
      p <- c(p, list(theme(legend.position='none')))
  }
  
  # Contour representation
  if (! is.null(name_contour))
  {
    p = c(p, gridContour(dbgrid, name=name_contour, useSel=useSel, ...))
    
    # Set the title                    
    title = paste(title, name_contour, sep=" ")
    
    # Set the Legend
    if (show.legend.contour)
    {
      if (is.null(legend.name.contour)) legend.name.contour = name_contour
      p <- c(p, list(labs(contour = legend.name.contour)))
    }
    else
      p <- c(p, list(theme(legend.position='none')))
  }  
  
  # Palette definition
  if (! is.null(palette))
    p <- c(p, .scaleColorFill(palette, na.value=na.value, limits=limits, ...))
  
  # Decoration
  p <- c(p, plot.decoration(title = title))
  
  p
}

#' Representing a Polygon
#' @param poly A Polygon object from the gstlearn library
#' @param show.title Plot the decoration attached to the figure
#' @param ... List of arguments passed to geom_polygon( )
#' @return The ggplot object
plot.polygon <- function(poly, show.title=FALSE, ...)
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
    p <- c(p, do.call(geom_polygon, c(list(data = df, mapping=aes(x=x,y=y)),  dotloc)))
  }
  
  # Decoration
  if (show.title)
    p <- c(p, plot.decoration(title = paste("Number of Polygons = ",npol)))
  
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
  
  p <- c(p, geom_histogram(data=df, mapping=aes(x=val), na.rm=TRUE, ...))
  p <- c(p, plot.decoration(title=name, xlab="Value", ylab="Count"))
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
  
  p <- c(p, geom_histogram(data = df, mapping=aes(x=val), na.rm=TRUE, ...))
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
  
  p <- c(p, geom_line(data = df, mapping=aes(x=absc,y=data), na.rm=TRUE, ...))
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
    p <- c(p, geom_path(data = df, mapping=aes(x=x,y=y), na.rm=TRUE, ...))
  
  if (flagPoint)
    p <- c(p, geom_point(data = df, mapping=aes(x=x,y=y), na.rm=TRUE, ...))
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
  p = c(p, geom_bin2d(data = df, mapping = aes(x=x, y=y), na.rm=TRUE, ...))
  p = c(p, list(theme_bw()))
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
  
  p = c(p, plot.XY(res$getY(), res$getZ(), ...))
  p = c(p, plot.geometry(xlim=res$getAylim(), ylim=res$getAzlim()))
  p = c(p, plot.decoration(xlab = "Gaussian", ylab = "Raw"))
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
#' @param diag_color Color of the diagonal
#' @param diag_line Line type of the diagonal
#' @param flagRegr Represent the linear regression of Y|X 
#' @param regr_color Color of the linear regression of Y|X
#' @param regr_line Line type of the linear regression of Y|X
#' @param flagBiss Represent the first bisector (Y=X) 
#' @param biss_color Color of the first bisector (Y=X)
#' @param biss_line Line type of the first bisector (Y=X)
#' @param flagSameAxes Define the same bounds for horizontal and vertical axes
#' @param show.legend.raster Show the legend when representing grid of occurrences (asPoint = FALSE)
#' @param legend.name.raster Name of the legend when representing grid of occurrences (asPoint = FALSE)
#' @param ... List of arguments passed to plot.XY() or plot.hist2d()
#' @return The ggplot object
plot.correlation <- function(db1, namex, namey, db2=NULL, useSel=TRUE,
    asPoint=FALSE, 
    flagDiag=FALSE, diag_color = "red", diag_line = "solid", 
    flagRegr=FALSE, regr_color = "blue", regr_line = "solid", 
    flagBiss=FALSE, biss_color = "green", biss_line = "solid", 
    flagSameAxes=FALSE, show.legend.raster = FALSE, legend.name.raster = NULL,
    ...)
{
  if (is.null(db2)) db2 = db1
  x = db1$getColumn(namex, useSel)
  y = db2$getColumn(namey, useSel)
  
  p = list()
  if (asPoint)
    p = c(p, plot.XY(x, y, flagLine=FALSE, flagPoint=TRUE, ...))
  else
    p = c(p, plot.hist2d(x, y, ...))
  
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
    p <- c(p, geom_segment(aes(x=xmin,y=ymin,xend=xmax,yend=ymax),
        linetype = diag_line, color = diag_color, na.rm=TRUE))
  }
  
  if (flagBiss)
  {
    bmin = min(xmin, ymin)
    bmax = max(xmax, ymax)
    p <- c(p, geom_segment(aes(x=bmin,y=bmin,xend=bmax,yend=bmax),
        linetype = biss_line, color = biss_color, na.rm=TRUE))
  }
  
  if (flagRegr)
  {
    regr = regression(db2, namey, namex, flagCste=TRUE)
    if (regr$nvar > 0)
    {
      a = regr$coeffs[1]
      b = regr$coeffs[2]
      ymin = a + b * xmin
      ymax = a + b * xmax
      p <- c(p, geom_segment(aes(x=xmin,y=ymin,xend=xmax,yend=ymax),
          linetype = regr_line, color = regr_color, na.rm=TRUE))
    }
  }
  
  p = c(p, plot.decoration(xlab=namex, ylab=namey))
  
  # Set the Legend
  if (show.legend.raster)
  {
    if (is.null(legend.name.raster)) legend.name.raster = "Count"
    p <- c(p, list(guides(fill = guide_colorbar(title=legend.name.raster, reverse=FALSE))))
  }
  else
    p <- c(p, list(theme(legend.position='none')))
  
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
  
  p = c(p, geom_rect(data = df, mapping=aes(xmin = xmin, xmax = xmax, 
          ymin = ymin, ymax = ymax, fill = colors), na.rm=TRUE, ...))
  p = c(p, plot.geometry(xlim=c(-maxG,+maxG), ylim=c(-maxG,+maxG)))
  p
}

#' Represent a meshing
#' @param mesh A mesh object from gstlearn library
#' @param flagFace Assign a color to each polygon of the meshing
#' @param flagApex Draw the rank of each apex of the meshing
#' @param ... List of arguments passed to geom_polygon(), geom_path() or geom_point()
#' @return The ggplot object
plot.mesh <- function(mesh, flagFace=FALSE, flagApex=FALSE, ...)
{
  p = list()
  
  nmesh = mesh$getNMeshes()
  for (imesh in 1:nmesh)
  {
    x = mesh$getCoordinatesPerMesh(imesh-1, 0, TRUE)
    y = mesh$getCoordinatesPerMesh(imesh-1, 1, TRUE)
    df = data.frame(x, y)
    if (flagFace)
      p <- c(p, geom_polygon(data = df, mapping=aes(x=x,y=y), ...))
    else
      p <- c(p, geom_path(data = df, mapping=aes(x=x,y=y), ...))
    if (flagApex)
      p <- c(p, geom_point(data = df, mapping=aes(x=x, y=y)))
  }  
  p
}

#' Represent the neighborhood
#' @param neigh A Neigh object from the gstlearn library
#' @param grid The target Data base from the gstlearn library
#' @param node Rank of the target 
#' @param flagCell Represent the target as the corresponding cell
#' @param flagZoom Zoom to the extension of the neighborhood
#' @return The ggplot object
plot.neigh <- function(neigh, grid, node=0, flagCell=FALSE, flagZoom=FALSE, ...)
{
	p = list()

    # Identify target location
    target = grid$getSampleCoordinates(node)
    
    # Represent the target location
    p = c(p, plot.XY(target[1], target[2], flagLine=FALSE, flagPoint=TRUE, ...))
    
    # Represent the edge of the target (if block)
    edges = grid$getCellEdges(node)
    p = c(p, plot.XY(edges[[1]], edges[[2]], ...))
    
    # Represent the Neighborhood Ellipsoid
    if (neigh$getType()$getValue() == ENeigh_MOVING()$getValue())
    {
    	edges = neigh$getEllipsoid(target)
    	p = c(p, plot.XY(edges[[1]], edges[[2]], ...))
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
       		p = c(p, plot.XY(cx, cy, flagLine=TRUE, flagPoint=FALSE, ...))
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
#' @return The ggplot object
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
                 arrow = arrow(length = unit(0.1, "cm")))
    
	p
}

setMethod("plot", signature(x="_p_AMesh"), function(x,y=missing,...)   plot.mesh(x,...))
setMethod("plot", signature(x="_p_DbGrid"), function(x,y="missing",...)  plot.grid(x,...))

setMethod("plot", signature(x="_p_Db"), function(x,y=missing,...) plot.point(x,...))
setMethod("plot", signature(x="_p_Polygons"), function(x,y=missing,...) plot.polygon(x,...))

setMethod("plot", signature(x="_p_Vario"), function(x,y=missing,...) plot.vario(x,...))
setMethod("plot", signature(x="_p_Model"), function(x,y="missing",...) plot.model(x,...))

setMethod("plot", signature(x="_p_Rule"), function(x,y="missing",...) plot.rule(x,...))
setMethod("plot", signature(x="_p_AAnam"), function(x,y="missing",...) plot.anam(x,...))
