################################################################################
#                                                                              #
#                            gstlearn R package                                #
#                                                                              #
# Copyright (c) (2023) MINES PARIS / ARMINES                                   #
# Authors: gstlearn Team                                                       #
# Website: https://github.com/gstlearn                                         #
# License: GPL v3                                                              #
#                                                                              #
################################################################################
#
# This is a set of functions working with ggplot2() which enable performing
# plots easily 
#

# Define the global values in the given environment position (search)
plot.initialize <- function(pos=1) 
{
  assign("plot.default_dims", list(c(8,8), c(8,8)), pos=pos)
  assign("plot.default_xlim", list(c(NA,NA), c(NA,NA)), pos=pos)
  assign("plot.default_ylim", list(c(NA,NA), c(NA,NA)), pos=pos)
  assign("plot.default_asp",  c(0, 1), pos=pos)
  invisible()
}

isNotDef <- function(arg)
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

plot.setDefaultGeographic <- function(dims=NA, xlim=NA, ylim=NA, asp=NA)
{
  plot.setDefaultInternal(2, dims=dims, xlim=xlim, ylim=ylim, asp=asp)
}

plot.setDefault <- function(dims=NA, xlim=NA, ylim=NA, asp=NA)
{
  plot.setDefaultInternal(1, dims=dims, xlim=xlim, ylim=ylim, asp=asp)
}

plot.setDefaultInternal <- function(mode=1, dims=NA, xlim=NA, ylim=NA, asp=NA)
{
  if (!isNotDef(dims))
    plot.default_dims[[mode]] = dims
  if (!isNotDef(xlim))
    plot.default_xlim[[mode]] = xlim
  if (!isNotDef(ylim))
    plot.default_ylim[[mode]] = ylim    
  if (!isNotDef(asp))
    plot.default_asp[[mode]] = asp
}

plot.printDefault <- function()
{
  for (mode in 1:2)
  {
    if (mode == 1)
      cat("Non geographical defaults (mode=1):\n")
    else
      cat("Geographical defaults (mode=2):\n")
    
    if (!isNotDef(plot.default_dims[[mode]]))
      cat("- Figure dimensions =", plot.default_dims[[mode]],"\n")
    else
      cat("- Figure dimensions (not defined)\n")
    
    if (!isNotDef(plot.default_xlim[[mode]]))
      cat("- Limits along X =",plot.default_xlim[[mode]],"\n")
    else
      cat("- Limits along X (not defined)\n")
    
    if (!isNotDef(plot.default_ylim[[mode]]))
      cat("- Limits along Y =",plot.default_ylim[[mode]],"\n")
    else
      cat("- Limits along Y (not defined)\n")
    
    if (plot.default_asp[mode] != 0)
      cat("- Aspect =",plot.default_asp[mode],"\n")
    else
      cat("- Aspect (automatic)\n")
  }    
}

scale_col_fill <- function(palette, ...)
{
  rcb <- c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", "PuBu",
      "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd",
      "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3",
      "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")
  rcb_num <- 1:18
  v <- c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo",
      "A", "B", "C", "D", "E", "F", "G", "H")
  
  if(length(palette) == 1) {
    if (any(palette == rcb) | any(palette == rcb_num)) {
      layer = scale_color_distiller(palette= palette, aesthetics= c("colour", "fill"), ...)
    } else if(any(palette == v)) {
     layer = scale_color_viridis_c(option= palette, aesthetics= c("colour", "fill"), ...)
    } 
  } else if(length(palette) == 2) {
    low = palette[1]
    high = palette[2]
    layer = scale_color_gradient(low= low, high= high, aesthetics= c("colour", "fill"), ...)
  } else if(length(palette) == 3) {
    low = palette[1]
    mid = palette[2]
    high = palette[3]
    layer = scale_color_gradient2(low= low, mid= mid, high= high, aesthetics= c("colour", "fill"), ...)
  } else {
    layer = scale_color_continuous(type= palette, aesthetics= c("colour", "fill"), ...)
  }
  layer
}

get.colors <- function()
{
  c("blue", "red", "green", "brown", "orange", "purple", "yellow")
}

#'@title ggPrint
#'@description Print the contents of a ggplot, possibly without warnings
#'@param p Current contents of the ggplot()
#'@param flag_suppress_warnings TRUE to suppress informational warnings
ggPrint <- function(p, flag_suppress_warnings = TRUE)
{
  if (flag_suppress_warnings)
    suppressWarnings(plot(p))
  else
    plot(p)
  invisible()
}

ggDefaultGeographic <- function()
{
  p <- ggplot()
  mode = 2
  p <- p + plot.geometry(dims=plot.default_dims[[mode]], 
        xlim=plot.default_xlim[[mode]], 
        ylim=plot.default_ylim[[mode]], 
        asp=plot.default_asp[mode])
  p
}

ggDefault <- function()
{
  p <- ggplot()
  mode = 1
  p <- p + plot.geometry(asp= plot.default_asp[mode])
  p
}

is_array <- function(arg, ndim=NA)
{
  if (length(arg) <= 1) return (FALSE)
  
  if (!isNotDef(ndim) && length(arg) != ndim) return (FALSE)
  
  TRUE
}

plot.decoration <- function(xlab = NA, ylab = NA, title = NA)
{
  p = list()
  if (!isNotDef(xlab))
    p <- append(p, list(labs(x = xlab)))
  if (!isNotDef(ylab))
    p <- append(p, list(labs(y = ylab)))
  if (!isNotDef(title))
  {
    p <- append(p, list(labs(title = title)))
    p <- append(p, list(theme(plot.title = element_text(hjust = 0.5))))
  }
  p
}

# Set the Geometry for the plot 'p'
# asp: Specify a value of "0" for an automatic aspect ratio
#
plot.geometry <- function(dims=NA, xlim=NA, ylim=NA, asp=NA, expand=waiver())
{
  p = list()
  if (! isNotDef(dims[1]))
  {
    if (is_array(dims, 2))
    {
      options(repr.p.width  = dims[1], repr.p.height = dims[2])
    }
    else
      cat("'dims' should be [a,b]. Ignored\n")
  }
  
  if (length(xlim) > 1)
  {
    if (is_array(xlim, 2))
    {
      p <- c(p, scale_x_continuous(limits=xlim, expand=expand))
    }
    else
      cat("'xlim' should be [a,b] or [NA,b] or [a,NA]. Ignored\n")
  }
  
  if (length(ylim) > 1)
  {
    if (is_array(ylim, 2))
    {
      p <- c(p, scale_y_continuous(limits=ylim, expand=expand))
    }
    else
      cat("'ylim' should be [a,b] or [NA,b] or [a,NA]. Ignored\n")
  }
  
  if (!isNotDef(asp))
  {     
    if (asp != 0)
      p = c(p, coord_fixed(asp))
    else
      p = c(p, list(theme(aspect.ratio = 1)))
  }
  
  p
}

# Function for representing a Model
plot.model <- function(model, ivar=0, jvar=0, vario=NA, idir=0, ...)
{
  p = list()
  p = c(p, plot.varmod(vario=vario, model=model, ivar=ivar, jvar=jvar, idir=idir,
          draw_vario=FALSE, ...))
  p = c(p, plot.decoration(xlab = "Distance", ylab = "Variogram"))
  p  
}

# Function for representing the Experimental Variogram together with the Model (optional)

plot.vario <- function(vario, ivar=0, jvar=0, idir=0,...)
{
  p = list()
  p = c(p, plot.varmod(vario=vario, ivar=ivar, jvar=jvar, idir=idir,...))
  p
}

selectItems <- function(nvalues, sitem=-1)
{
  if (sitem >= 0)
    outs = sitem
  else
    outs = seq(0, nvalues-1)
  outs
}

varioElem <- function(vario, ivar=0, jvar=0, idir=0, 
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

modelElem <- function(model, ivar=0, jvar=0, codir=NA,
    nh = 100, hmax = NA, asCov=FALSE, flag_envelop = TRUE, 
    env_color='black', env_linetype="dashed", env_size=0.5,
    ...)
{
  dots = list(...)
  p = list()
  ndim = model$getDimensionNumber()
  
  # if hmax not specified = 3*maximum range of the model's basic structures
  if (isNotDef(hmax))
  {
    hmax = 0
    for (icova in 1:model$getCovaNumber())
    {
      range_max = max(model$getCova(icova-1)$getRanges())
      if (3*range_max > hmax)
        hmax = 3*range_max
    }
  }
  if (isNotDef(hmax) || hmax == 0) hmax = 1.
  
  if (isNotDef(codir))
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
  gg = model$sample(hh, ivar=ivar, jvar=jvar, codir=codir,
      nostd=0, asCov=asCov)
  df = data.frame(gg = gg[istart:nh], hh = hh[istart:nh])
  p = c(p, geom_line(data = df, mapping=aes(x=hh, y=gg), na.rm=TRUE, ...))
  
  # Represent the coregionalization envelop
  if (ivar != jvar && flag_envelop)
  {
    gg = model$sample(hh, ivar=ivar, jvar=jvar, codir=codir, 
        nostd=-1, asCov=asCov)
    df = data.frame(gg = gg[istart:nh], hh = hh[istart:nh])
    p = c(p, geom_line(data = df, mapping=aes(x=hh, y=gg), na.rm=TRUE, 
            color = env_color, linetype = env_linetype, size=env_size))
    
    gg = model$sample(hh, ivar=ivar, jvar=jvar, codir=codir, 
        nostd=1, asCov=asCov)
    df = data.frame(gg = gg[istart:nh], hh = hh[istart:nh])
    p = c(p, geom_line(data = df, mapping=aes(x=hh, y=gg), na.rm=TRUE, 
            color = env_color, linetype = env_linetype, size=env_size))
  }
  p
}

# This function must not be added to a prior call to ggplot().
#
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
  if (! isNotDef(vario)) ndir = vario$getDirectionNumber()
  nvar = 1
  if (! isNotDef(vario)) nvar = vario$getVariableNumber()
  if (! isNotDef(model)) nvar = model$getVariableNumber()
  cols = get.colors()
  
  idirUtil = selectItems(ndir, idir)
  ivarUtil = selectItems(nvar, ivar)
  jvarUtil = selectItems(nvar, jvar)
  ivarN = length(ivarUtil)
  jvarN = length(jvarUtil)
  
  if (isNotDef(hmax))
  {
    if (! isNotDef(vario))
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
        if (! isNotDef(vario) && draw_vario)
        {
          dotloc = dots
          if (! has_color) dotloc$color=cols[idir+1]
          if (! has_linetype) dotloc$linetype = "dashed"
          p = c(p, do.call(varioElem, c(list(vario=vario, ivar=ivar, jvar=jvar, idir=idir, 
                          var_color=var_color, var_linetype=var_linetype, var_size=var_size,
                          draw_variance=draw_variance, draw_psize=draw_psize, 
                          draw_plabel=draw_plabel, label=label, flagLimits=FALSE), dotloc)))
        }
        
        # Plotting the Model
        if (! isNotDef(model))
        {
          dotloc = dots
          if (! has_color) dotloc$color=cols[idir+1]
          if (! has_linetype) dotloc$linetype = "solid"
          if (! isNotDef(vario) && ! has_codir) dotloc$codir = vario$getCodirs(idir) 
          p = c(p, do.call(modelElem, c(list(model, ivar, jvar,  
                          nh = nh, hmax = hmax, asCov=asCov, flag_envelop=flag_envelop,
                          env_color = env_color, env_linetype = env_linetype, 
                          env_size=env_size), dotloc)))
        }
        
        # Adding some decoration
        p = c(p, plot.decoration(xlab = "Distance", ylab = "Variogram"))
      }
      
      # Informing bound criterion
      if (! isNotDef(vario))
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

# Arrange a set of figures for the multivariate representation
multi.varmod <- function(vario, model=NA, ivar=-1, jvar=-1, idir=-1,
    nh = 100, hmax = NA, draw_psize=-1, draw_plabel=FALSE, 
    asCov=FALSE, draw_variance = TRUE, flag_envelop=TRUE, 
    var_color='black', var_linetype="dashed", var_size=0.5, 
    env_color='black', env_linetype="dashed", env_size=0.5,
    label=NULL, ...)
{
  nvar = vario$getVariableNumber()
  
  ivarUtil = selectItems(nvar, ivar)
  jvarUtil = selectItems(nvar, jvar)
  ivarN = length(ivarUtil)
  jvarN = length(jvarUtil)
  
  if (isNotDef(hmax))
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

readPointCoor <- function(db)
{
  x = db$getCoordinates(0,TRUE)
  y = db$getCoordinates(1,TRUE)
  df = data.frame(x,y)
  df
}

readGridCoor <- function(dbgrid, name, usesel= FALSE)
{
  x = dbgrid$getColumnByLocator(ELoc_X(),0, FALSE, FALSE)
  y = dbgrid$getColumnByLocator(ELoc_X(),1, FALSE, FALSE)
  data = dbgrid$getColumn(name, usesel, FALSE)
  if (length(data) != length(x))
  {
    cat("Variable",name,"does not exist or does not have correction dimension\n")
    stop()
  }
  df = data.frame(x,y,data)
  df
}

# Function for plotting a point data base, with optional color and size variables
pointSymbol <- function(db, name_color=NULL, name_size=NULL,
    flagAbsSize = FALSE, flagCst=FALSE,
    ...) 
{ 
  # Creating the necessary data frame
  df = readPointCoor(db)
  
  # Color of symbol
  colval = NULL
  if (! is.null(name_color)) {
    colval  = db$getColumn(name_color, TRUE)
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

# Function for plotting a point data base, with label variables
pointLabel <- function(db, name, digit=2, ...) 
{  
  # Creating the necessary data frame
  df = readPointCoor(db)
  
  # Label of symbols
  labval  = round(db$getColumn(name,TRUE),digit)
  df["labval"] = as.character(labval)
  
  layer <- geom_text(data = df, mapping=aes(x=x, y=y, label=labval), 
      na.rm=TRUE, ...)
  
  layer
}

# Function for plotting a point data base, with optional color and size variables
plot.point <- function(db, name_color=NULL, name_size=NULL, name_label=NULL,
    sizmin=1, sizmax=5, flagAbsSize = FALSE, flagCst=FALSE,
    show.legend.symbol=FALSE, legend.name.color="P-Color",
    legend.name.size="P-Size",
    show.legend.label=FALSE, legend.name.label="P-Label", ...)
{ 
  p = list()
  title = ""
  
# If no variable is defined, use the default variable for Symbol(size) representation
# The default variable is the first Z-locator one, or the last variable in the file
  flagTitleDefault = FALSE
  if (is.null(name_color) && is.null(name_size) && is.null(name_label))
  {
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
        flagAbsSize = flagAbsSize, flagCst=flagCst,
        show.legend = show.legend.symbol,
        ...))
    
    if (! is.null(name_size) && ! flagCst)
      p <- c(p, scale_size_continuous(range = c(sizmin, sizmax)))
    
    # Set the default title
    if (! is.null(name_color))
      title = paste(title, name_color)
    if (! is.null(name_size))
      title = paste(title, name_size, sep=" ")
    
    # Set the Legend
    p <- c(p, list(labs(color = legend.name.color)))
    p <- c(p, list(labs(size = legend.name.size)))
  }
  
  if (! is.null(name_label))
  {
    p <- c(p, pointLabel(db, name=name_label, digit=digit_label, 
        show.legend=show.legend.label, ...))
    
    # Set the title              
    title = paste(title, name_label, sep=" ")
    
    # Set the legend
    p <- c(p, list(labs(label = legend.name.label)))
  }
  
  # Decoration
  if (flagTitleDefault) title = "Sample Location"
  p <- c(p, plot.decoration(title = title))
  
  p
}

gridRaster <- function(dbgrid, name, usesel = TRUE, ...)
{
  # Reading the Grid information
  df = readGridCoor(dbgrid, name, usesel)
  
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

gridContour <- function(dbgrid, name, usesel = TRUE, ...)
{
  # Reading the Grid information
  df = readGridCoor(dbgrid, name, usesel)
  
  layer <- geom_contour(data = df, mapping=aes(x = x, y = y, z = data), ...)
  
  layer
}

#
# Function for plotting a variable informed in a grid Db
#
plot.grid <- function(dbgrid, name_raster=NULL, name_contour=NULL,
    usesel = TRUE, palette=NULL, na.value = "white", limits = NULL, 
    show.legend.raster=FALSE, legend.name.raster="G-Raster", 
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
    if (dbgrid$getLocNumber(ELoc_Z()) > 0)
      name_raster = dbgrid$getNameByLocator(ELoc_Z(),0)
    else
      # if no Z locator, choose the last field
      name_raster = dbgrid$getLastName()
  }
  
  # Raster representation
  
  if (! is.null(name_raster))
  {
    p <- c(p, gridRaster(dbgrid, name=name_raster, usesel=usesel, ...))
    
    # Set the title
    title = paste(title,name_raster)
    
    # Set the Legend
    if (show.legend.raster)
      p <- c(p, list(guides(fill = guide_colorbar(title=legend.name.raster, reverse=FALSE))))
    else
      p <- c(p, list(theme(legend.position='none')))
  }
  
  # Contour representation
  
  if (! is.null(name_contour))
  {
    p = c(p, gridContour(dbgrid, name=name_contour, usesel=usesel, ...))
    
    # Set the title                    
    title = paste(title, name_contour, sep=" ")
  }  
  
  # Palette definition

  if (! is.null(palette))
    p <- c(p, scale_col_fill(palette, na.value=na.value, limits=limits, ...))
  
  # Decoration
  p <- c(p, plot.decoration(title = title))
  
  p
}

# Function to display a polygon (not tested)

plot.polygon <- function(poly, show.title=FALSE, ...)
{
  dots = list(...)
  has_color = "color" %in% names(dots)
  
  p = list()
  npol = poly$getPolySetNumber()
  cols = get.colors()
  
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

# Function for plotting the histogram of a variable
plot.hist <- function(db, name, usesel=TRUE, ...)
{
  p = list()
  val  = db$getColumn(name, usesel)
  df = data.frame(val)
  
  p <- c(p, geom_histogram(data=df, mapping=aes(x=val), na.rm=TRUE, ...))
  p <- c(p, plot.decoration(title=name, xlab="Value", ylab="Count"))
  p
}

# Function for plotting histogram for a table of values
plot.histArray <- function(val, ...)
{
  p = list()  
  df = data.frame(val)
  
  p <- c(p, geom_histogram(data = df, mapping=aes(x=val), na.rm=TRUE, ...))
  p
}

# Function for plotting a curve of regularly sampled values
plot.curve <- function(data, ...)
{
  p = list()  
  absc = seq(1,length(data))
  df = data.frame(absc,data)
  
  p <- c(p, geom_line(data = df, mapping=aes(x=absc,y=data), na.rm=TRUE, ...))
  p
}

# Function for representing a set of points with optional symbols and joining lines
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

# Function for representing an anamorphosis
plot.anam <- function(anam, ndisc=100, aymin=-10, aymax=10, ...)
{
  p = list()
  res = anam$sample(ndisc, aymin, aymax)
  
  p = c(p, plot.XY(res$getY(), res$getZ(), ...))
  p = c(p, plot.geometry(xlim=res$getAylim(), ylim=res$getAzlim()))
  p = c(p, plot.decoration(xlab = "Gaussian", ylab = "Raw"))
  p
}

# Function for representing a scatter plot
plot.correlation <- function(db1, name1, name2, db2=NULL, usesel=FALSE,
    asPoint=FALSE, 
    flagDiag=FALSE, diag_color = "red", diag_line = "solid", 
    flagRegr=FALSE, regr_color = "blue", regr_line = "solid", 
    flagBiss=FALSE, biss_color = "green", biss_line = "solid", 
    flagSameAxes=FALSE, 
    show.legend.raster = FALSE, legend.name.raster="Count",
    ...)
{
  if (is.null(db2)) db2 = db1
  x = db1$getColumn(name1, usesel)
  y = db2$getColumn(name2, usesel)
  
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
    regr = regression(db2, name2, name1, flagCste=TRUE)
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
  
  p = c(p, plot.decoration(xlab=name1, ylab=name2))
  
  # Set the Legend
  if (show.legend.raster)
    p <- c(p, list(guides(fill = guide_colorbar(title=legend.name.raster, reverse=FALSE))))
  else
    p <- c(p, list(theme(legend.position='none')))
  
  p 
}

# Representing a Lithotype rule
plot.rule <- function(rule, proportions=NULL, maxG = 3., ...)
{
  p = list()
  nrect = rule$getFaciesNumber()
  if (! is.null(proportions)) 
    rule$setProportions(proportions)
  else
    rule$setProportions()
  cols = get.colors()
  
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


# Function to display a polygon (not tested)
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

setMethod("plot", signature(x="_p_AMesh"), function(x,y=missing,...)   plot.mesh(x,...))
setMethod("plot", signature(x="_p_DbGrid"), function(x,y="missing",...)  plot.grid(x,...))

setMethod("plot", signature(x="_p_Db"), function(x,y=missing,...) plot.point(x,...))
setMethod("plot", signature(x="_p_Polygons"), function(x,y=missing,...) plot.polygon(x,...))

setMethod("plot", signature(x="_p_Vario"), function(x,y=missing,...) plot.vario(x,...))
setMethod("plot", signature(x="_p_Model"), function(x,y="missing",...) plot.model(x,...))

setMethod("plot", signature(x="_p_Rule"), function(x,y="missing",...) plot.rule(x,...))
setMethod("plot", signature(x="_p_AAnam"), function(x,y="missing",...) plot.anam(x,...))
