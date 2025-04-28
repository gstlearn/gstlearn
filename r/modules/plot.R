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

#' Check if an argument is defined
#' @param arg Argument to be checked
#' 
#' @keywords internal
#' @noRd
.isNotDef <- function(arg) {
  if (is.null(arg)) {
    return(TRUE)
  }
  warn.old = options("warn")
  options(warn = -1)
  if (length(arg) == 1) {
    if (is.na(arg)) {
      options(warn.old)
      return(TRUE)
    }
  } else {
    for (i in 1:length(arg))
    {
      if (is.na(arg[i])) {
        options(warn.old)
        return(TRUE)
      }
    }
  }
  options(warn.old)
  return(FALSE)
}

#' Function to create breaks (used to plot size legend
#' 
#' @keywords internal
#' @noRd
.getSeq <- function(v, n = 4, nd = 3) {
  u = floor(seq(from = floor(v[1] * 10**nd), to = ceiling(v[2] * 10**nd), length.out = n)) / 10**nd
  u = pretty(v, n)
  return(u)
}

#' Initialize a new figure and define its environment
#'
#' Note: This method includes the initial call to ggplot()
#'
#' @param dims Vector giving the dimensions of the figure
#' @param xlim Bounds of the figure along the horizontal axis (when left to NA, it will be adjusted to the figure contents)
#' @param ylim Bounds of the figure along the vertical axis (when left to NA, it will be adjusted to the figure contents)
#' @param asp Aspect ratio Y/X. A value of 1 means that the scales will be similar on the two axes of the figure.
#'
plot.init <- function(dims = NA, xlim = NA, ylim = NA, asp = NA) {
  if (!require(ggplot2, quietly = TRUE)) {
    stop("Package 'ggplot2' is mandatory to use this function!")
  }

  # Define the size of the graphics
  if (.isArray(dims, 2)) {
    options(repr.plot.width = dims[1], repr.plot.height = dims[2])
  }

  # Initialize a new plot
  p = ggplot()

  # Set the geometry of the figure
  p = p + plot.geometry(xlim = xlim, ylim = ylim, asp = asp)

  # Return the ggplot object
  p
}

#' Allow redefining a new aesthetic element (if already defined)
#' @param aestype Should be "colour" or "fill" or "linetype" or "size"
#' 
#' @keywords internal
#' @noRd
.appendNewScale <- function(p, aestype) {
  p <- append(p, list(new_scale(aestype)))
  p
}

#' Generate the list of predefined color palettes
#'
#' @return A data Frame containing the names and origins or the different palettes
#'
#' @keywords internal
#' @noRd
.getAllPalettes <- function() {
  if (!require(RColorBrewer, quietly = TRUE)) {
    stop("Package 'RColorBrewer' is mandatory to use this function!")
  }

  rcb <- rownames(brewer.pal.info)
  rcb_num <- 1:18
  v <- c(
    "magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo",
    "A", "B", "C", "D", "E", "F", "G", "H"
  )

  mat = rbind(
    cbind(rcb, "brewer"),
    cbind(v, "viridis"),
    cbind(rcb_num, "brewer")
  )
  df = data.frame(name = mat[, 1], origin = mat[, 2])
  df
}

#' Define the "colour" using input 'palette' definition
#' @param palette Reference palette used for defining the current color map
#' @param naColor Color used for representing NA values
#' @param flagDiscrete True for defining a Discrete Color scale
#' @param limits Limits for the color scale
#' @param title Title of the color scale
#' 
#' @keywords internal
#' @noRd
.defineColour <- function(palette, naColor = "transparent", flagDiscrete = FALSE, limits = NULL, title = NA) {
  dfPalettes = .getAllPalettes()
  aes_list = c("color")
  name = .defineName(title)

  if (flagDiscrete) {
    # Case of a Discrete color scale
    layer = scale_colour_manual(
      values = palette, aesthetics = aes_list,
      na.value = naColor, name = name
    )
  } else {
    # Case of a Continuous Color Scale
    if (length(palette) == 0) {
      # Case where 'palette' is not provided
      layer = scale_colour_viridis_c(
        option = "viridis", aesthetics = aes_list,
        na.value = naColor, limits = limits, name = name
      )
    } else if (length(palette) == 1) {
      # Palette is provided: find its rank in the list of recorded palettes
      rank = which(dfPalettes$name == palette)

      if (length(rank) == 0) {
        # Case where 'palette' is not present in the palette list
        cat("Palette (",palette,") is not defined in the list of palettes: 'viridis' is used instead")
        layer = scale_colour_viridis_c(
          option = "viridis", aesthetics = aes_list,
          na.value = naColor, limits = limits, name = name
        )
      } else {
        origin = dfPalettes$origin[rank]
        if (origin == "brewer") {
          # Case where 'palette' is a RColorBrewer palette
          layer = scale_colour_distiller(
            palette = palette, aesthetics = aes_list,
            na.value = naColor, name = name
          )
        } else if (origin == "viridis") {
          # Case where 'palette' is a viridis palette
          layer = scale_colour_viridis_c(
            option = palette, aesthetics = aes_list,
            na.value = naColor, limits = limits, name = name
          )
        } else {
          cat("Palette origin does not exist")
        }
      }
    } else {
      layer = scale_colour_gradientn(
        colours = palette, aesthetics = aes_list,
        na.value = naColor, limits = limits, name = name
      )
    }
  }
  return(layer)
}

#' Define the "fill" using input 'palette' definition
#' @param palette Reference palette used for defining the current color map
#' @param naColor Color used for representing NA values
#' @param flagDiscrete True for defining a Discrete Color scale
#' @param limits Limits for the color scale
#' @param title Title of the color scale
#' 
#' @keywords internal
#' @noRd
.defineFill <- function(palette = NULL, naColor = "transparent", flagDiscrete = FALSE, limits = NULL, title = NA) {
  dfPalettes = .getAllPalettes()
  aes_list = c("fill")
  name = .defineName(title)

  if (flagDiscrete) {
    # Case of a Discrete color scale
    layer = scale_fill_manual(values = palette, na.value = naColor, name = name)
  } else {
    # Case of a Continuous Color Scale
    if (length(palette) == 0) {
      # Case where 'palette' is not provided
      layer = scale_fill_viridis_c(
        option = "viridis", aesthetics = aes_list,
        na.value = naColor, limits = limits, name = name
      )
    } else if (length(palette) == 1) {
      # Palette is provided: find its rank in the list of recorded palettes
      rank = which(dfPalettes$name == palette)

      if (length(rank) == 0) {
        # Case where 'palette' is not present in the palette list
        layer = scale_fill_gradientn(
          colours = palette, aesthetics = aes_list,
          na.value = naColor, limits = limits, name = name
        )
      } else {
        origin = dfPalettes$origin[rank]
        if (origin == "brewer") {
          layer = scale_fill_distiller(
            palette = palette, aesthetics = aes_list,
            na.value = naColor, limits = limits, name = name
          )
        } else if (origin == "viridis") {
          layer = scale_fill_viridis_c(
            option = palette, aesthetics = aes_list,
            na.value = naColor, limits = limits, name = name
          )
        } else {
          cat("Palette origin does not exist")
        }
      }
    } else {
      layer = scale_fill_gradientn(
        colours = palette, aesthetics = aes_list,
        na.value = naColor, limits = limits, name = name
      )
    }
  }
  return(layer)
}

#' Define a name by default
#' 
#' @keywords internal
#' @noRd
.defineName <- function(title = NA) {
  name = waiver()
  if (!is.na(title)) {
    if (is.character(title)) {
      name = title
    }
  }
  return(name)
}

#' Define a new scale
#' 
#' @keywords internal
#' @noRd
.defineSize <- function(sizval, sizeRange, title=NA)
{
  name = .defineName(title)
  breaks = .getSeq(range(sizval,na.rm=TRUE))
  layer = scale_size(breaks = breaks,limits=range(breaks),range = sizeRange, name=name)

  return(layer)
}

#' Define a series of distinct colors
#' 
#' @keywords internal
#' @noRd
.getColors <- function()
{
	c("blue", "red", "green", "brown", "orange", "purple", "yellow")
}


#' Check if the argument can be considered as an array (with possibly imposed dimensions)
#' @param arg Input argument
#' @param ndim Imposed dimension for the input argument (no check is performed if NA)
#' 
#' @keywords internal
#' @noRd
.isArray <- function(arg, ndim=NA)
{
  if (length(arg) <= 1) return (FALSE)
  if (.isNotDef(arg)) return (FALSE)
  if (length(arg) != ndim) return (FALSE)
  
  TRUE
}

#' Define the variable used by default (when not explictly defined)
#' @param db Data base from gstlearn
#' @param name1 Proposed input name
#' @param name2 Auxiliary input variable name
#'
#' @note The defaulted variable is the first one attached to the locator-Z (if any);
#' @note otherwise it is the last defined variable within 'db'.
#' @return Name of the defaulted variable
#' 
#' @keywords internal
#' @noRd
.defaultVariable <- function(db, name1 = NULL, name2 = NULL)
{
  # If at least one variable name is defined, return 'name1 as is'
  if (!is.null(name1) || !is.null(name2))
    return(name1)

  # If no variable is defined, return a variable name by default
  if (db$getNLoc(ELoc_Z()) > 0)
    # Use the first Z-Locator variable (if defined)
		return(db$getNameByLocator(ELoc_Z(),0))
  else
    # Use the last defined variable
	  return(db$getLastName())

  # No solution is found, return an error
  cat("No variable name is provided. Procedure is aborted")
  stop()
}

#' Check if the Legend must be defined or not
#' 
#' @param flagLegend Flag for displaying the legend
#' 
#' @keywords internal
#' @noRd
.showLegend <- function(flagLegend) {
  show = ifelse(!flagLegend, FALSE, NA)
  show
}

#' Select a list of items
#' @param nvalues Number of items in the list
#' @param sitem   Default value for the item withon this list (if non-negative)
#' @return The returned list of items 
#' 
#' @keywords internal
#' @noRd
.selectItemsInList <- function(nvalues, sitem=-1)
{
  if (sitem >= 0)
    outs = sitem
  else
    outs = seq(0, nvalues-1)
  outs
}

#' Draw an elementary experimental variogram
#' 
#' @keywords internal
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

#' Draw an elementary Model
#' 
#' @keywords internal
#' @noRd
.modelElementary <- function(model, ivar=0, jvar=0, codir=NA,
    nh = 100, hmax = NA, asCov=FALSE, flagEnvelop = TRUE, 
    envColor='black', envLinetype="dashed", envSize=0.5, 
    dirName = "Direction1", lineName="Model",
    ...)
{
  dots = list(...)
  p = list()
  ndim = model$getNDim()
  
  # if hmax not specified = 3*maximum range of the model's basic structures
  if (.isNotDef(hmax))
  {
    hmax = 0
    for (icova in 1:model$getNCov())
    {
      range_max = max(model$getCovAniso(icova-1)$getRanges())
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
  for (icova in 1:model$getNCov())
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
  gg = model$sample(hh, codir=codir, ivar=ivar, jvar=jvar, mode=mode)
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

  # Adding the vertical axis at X=0
  p = append(p, geom_vline(xintercept = 0., color='black', size=0.5))
  
  # Adding the horizontal axis at Y=0
  p = append(p, geom_hline(yintercept = 0., color = "black", size = 0.5))
  
  p
}

#' Read a set of sample coordinates
#' @param db Db item from the gstlearn library
#' @param useSel Use of an optional selection (masking off samples)
#' @param posX Rank of the coordinate which will serve as first coordinate
#' @param posY Rank of the coordinate which will serve as the second coordinate
#' @return a Dataframe containing the 2-D coordinates
#' 
#' @keywords internal
#' @noRd
.readPointCoor <- function(db, useSel=TRUE, posX=0, posY=1)
{
  if (db$getNDim() > 0) 
    x = db$getOneCoordinate(posX,useSel)
  if (db$getNDim() > 1)
    y = db$getOneCoordinate(posY,useSel)
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
#'
#' @keywords internal
#' @noRd
.readGridCoor <- function(dbgrid, name, useSel = FALSE, posX = 0, posY = 1, corner = NA) {
  if (.isNotDef(corner)) {
    corner = rep(0, dbgrid$getNDim())
  }

  if (dbgrid$getNDim() == 1) {
    data = dbgrid$getColumn(name, useSel, FALSE)
    x = dbgrid$getColumnByLocator(ELoc_X(), posX, FALSE, FALSE)
    y = dbgrid$getColumnByLocator(ELoc_X(), posY, FALSE, FALSE)
  } else {
    data = dbgrid$getOneSlice(name, posX, posY, corner, useSel)
    nameX = dbgrid$getNameByLocator(ELoc_X(), posX)
    x = dbgrid$getOneSlice(nameX, posX, posY, corner, FALSE)
    nameY = dbgrid$getNameByLocator(ELoc_X(), posY)
    y = dbgrid$getOneSlice(nameY, posX, posY, corner, FALSE)
  }

  if (length(data) != length(x)) {
    cat("Variable", name, "does not exist or does not have correction dimension\n")
    stop()
  }
  df = data.frame(x, y, data)
  df
}

#' Check if the input 'dbgrid' is a DbGrid.
#' Also check for rotation (optional)
#' 
#' @keywords internal
#' @noRd
.isGrid <- function(dbgrid, flagNoRotate = FALSE) 
{
  if (!dbgrid$isGrid()) {
    cat("This function is restricted to Grid Db and cannot be used here\n")
    return(FALSE)
  }
  if (flagNoRotate) 
  {
    if (dbgrid$isGridRotated()) 
    {
      cat("This function cannot be used with rotated grids\n")
      return(FALSE)
    }
  }
  return(TRUE)
}

#' Print the list of predefined color palettes in \pkg{gstlearn}.
#'
#' @return Prints the list of color palette names and returns nothing.
#'
printAllPalettes <- function() {
  print(.getAllPalettes())
}

#' Print the contents of a ggplot, possibly without warnings
#'
#' @param p Current contents of the ggplot()
#' @param flagSuppressWarnings TRUE to suppress informational warnings
plot.end <- function(p, flagSuppressWarnings = TRUE)
{
  if (flagSuppressWarnings)
    suppressWarnings(plot(p))
  else
    plot(p)
  invisible()
}

#' Draw the decoration of a figure (title, axis labels, ...)
#'
#' @param title Title of the figure
#' 
#' @param xlab Label along the horizontal axis
#' @param ylab Label along the vertical axis
#' @param flagDefaultTitle TRUE if argument 'title' must be set only if no title already defined
#' @return The ggplot object
plot.decoration <- function(title=NA, xlab = NA, ylab = NA, flagDefaultTitle = FALSE)
{
  p = list()
  if (!.isNotDef(xlab))
    p <- append(p, list(labs(x = xlab)))
  if (!.isNotDef(ylab))
    p <- append(p, list(labs(y = ylab)))
  if (!.isNotDef(title))
  {
    flagPlotTitle = TRUE
    if (flagDefaultTitle && ! is.null(p$labels$title)) flagPlotTitle = FALSE
    if (flagPlotTitle)
    {
        p <- append(p, list(labs(title = title)))
    	p <- append(p, list(theme(plot.title = element_text(hjust = 0.5))))
    }
  }
  p
}

#' Set the Geometry for the current plot
#'
#' @param xlim Bounds along the horizontal axis
#' @param ylim Bounds along the vertical axis
#' @param asp  Aspect Ratio Y/X ("0" for an automatic aspect ratio)
#' @param expand Adding padding around data
#' @return The ggplot object
plot.geometry <- function(xlim = NA, ylim = NA, asp = NA, expand = waiver()) {
  p = list()

  if (.isArray(xlim, 2)) {
    p <- append(p, scale_x_continuous(limits = xlim, expand = expand))
  }

  if (.isArray(ylim, 2)) {
    p <- append(p, scale_y_continuous(limits = ylim, expand = expand))
  }

  if (!.isNotDef(asp)) {
    if (asp != 0) {
      p = append(p, coord_fixed(asp))
    } else {
      p = append(p, list(theme(aspect.ratio = 1)))
    }
  }
  p
}

#' Function for representing a Model
#' @param model An object of class Model from gstlearn
#' @param ivar Rank of the variable to be represented (-1 for all variables)
#' @param jvar Rank of the second variable to be represented in multivariate case (-1 for all variables)
#' @param vario An object of class Vario of gstlearn (optional)
#' @param idir Rank of the direction (-1 for all directions)
#' @param ... Arguments passed to plot.varmod()
#'
#' If 'vario' is defined, the calculation direction is given by the definition of direction 'idir' within 'vario'
#'
#' @return The ggplot object
plot.model <- function(model, ivar = 0, jvar = 0, vario = NA, idir = 0, ...) {
  p = list()
  p = append(p, plot.varmod(
    vario = vario, model = model, ivar = ivar, jvar = jvar, idir = idir,
    drawVario = FALSE, ...
  ))
  p = append(p, plot.decoration(xlab = "Distance", ylab = "Variogram"))
  p
}

#' Function for representing the Experimental Variogram
#' @param vario An object of the class Vario of gstlearn
#' @param ivar Rank of the variable to be represented (-1 for all variables)
#' @param jvar Rank of the second variable to be represented in multivariate case (_& for all variables)
#' @param idir Rank of the direction to be represented (-1 for all directions)
#' @param ... Arguments passed to plot.varmod()
#' @return The ggplot object
plot.vario <- function(vario, ivar=-1, jvar=-1, idir=-1,...)
{
  p = list()
  p = append(p, plot.varmod(vario=vario, ivar=ivar, jvar=jvar, idir=idir,...))
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
#' @param cols List of colors (optional)
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
    cols=NA, drawVario=TRUE, flagLegend=FALSE, ...)
{
  if (!require(ggnewscale, quietly=TRUE))
    stop("Package ggnewscale is mandatory to use this function!")

  dots = list(...)
  has_codir = "codir" %in% names(dots)
  linetypes = NULL
  if (! .isNotDef(vario) && drawVario) linetypes = c(linetypes, varioLinetype)
  if (! .isNotDef(model)) linetypes = c(linetypes, modelLinetype)
  
  p = list()
  ndir = 1
  if (! .isNotDef(vario)) ndir = vario$getNDir()
  nvar = 1
  if (! .isNotDef(vario)) nvar = vario$getNVar()
  if (! .isNotDef(model)) nvar = model$getNVar()
  if (missing(cols)) cols = .getColors()
  
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

  p <- append(p, scale_linetype_manual(name="Types", values = linetypes))
  p <- .appendNewScale(p, "linetype")
  p <- append(p, scale_color_manual(name = "Directions", values = cols))
  p <- .appendNewScale(p, "colour")
  
  # Loop on the variables
  flag_allow_negative_X = FALSE
  flag_allow_negative_Y = FALSE
  
  for (ivar in ivarUtil)
  {
    for (jvar in jvarUtil)
    {
      
      # Define the current plot
      for (idir in idirUtil)
      {
         if (! .isNotDef(vario))
         {
          angles = GeometryHelper_rotationGetAngles(vario$getCodirs(idir),TRUE)
          dirName = paste("Vario dir =", paste(round(angles,3), collapse=' '))
         }
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
#' @param ... Arguments passed to plot.varmod()
#' @return The ggplot object
multi.varmod <- function(vario, model=NA, ivar=-1, jvar=-1, idir=-1,
    nh = 100, hmax = NA, drawPsize=-1, drawPlabel=FALSE, 
    asCov=FALSE, drawVariance = TRUE, flagEnvelop=TRUE, 
    varioLinetype = "dashed", modelLinetype = "solid",
    varColor='black', varLinetype="dashed", varSize=0.5, 
    envColor='black', envLinetype="dashed", envSize=0.5,
    ...)
{
  if (!require(ggpubr, quietly = TRUE)) {
    stop("Package ggpubr is mandatory to use this function!")
  }

  nvar = vario$getNVar()

  ivarUtil = .selectItemsInList(nvar, ivar)
  jvarUtil = .selectItemsInList(nvar, jvar)
  ivarN = length(ivarUtil)
  jvarN = length(jvarUtil)

  if (.isNotDef(hmax)) {
    hmax = vario$getHmax(ivar, jvar, idir)
  }

  # Loop on the variables

  index = 0
  plot_lst <- vector("list", length = ivarN * jvarN)

  for (ivar in ivarUtil) {
    for (jvar in jvarUtil)
    {
      # Define the current plot
      index = index + 1

      g = plot.init()

      if (ivar < jvar) {
        g = g + theme(
          panel.background = element_rect(fill = "transparent"), # transparent panel bg
          plot.background = element_rect(fill = "transparent", color = NA), # transparent plot bg
          panel.grid.major = element_blank(), # remove major gridlines
          panel.grid.minor = element_blank(), # remove minor gridlines
          legend.background = element_rect(fill = "transparent"), # transparent legend bg
          legend.box.background = element_rect(fill = "transparent")
        ) # transparent legend panel
      } else {
        g = g + plot.varmod(
          vario = vario, model = model, ivar = ivar, jvar = jvar, idir = -1,
          nh = nh, hmax = hmax, drawPsize = drawPsize, drawPlabel = drawPlabel,
          asCov = asCov, drawVariance = drawVariance, flagEnvelop = flagEnvelop,
          varioLinetype = varioLinetype, modelLinetype = modelLinetype,
          varColor = varColor, varLinetype = varLinetype, varSize = varSize,
          envColor = envColor, envLinetype = envLinetype, envSize = envSize,
          ...
        )
      }
      plot_lst[[index]] <- g
    }
  }
  p = ggarrange(plotlist = plot_lst, nrow = ivarN, ncol = jvarN)
  p
}

#' Plotting a point data base where samples are displayed with different color and/or size
#' @param db Data Base containing the information to be displayed
#' @param nameColor Name of the variable to be represented in color
#' @param nameSize Name of the variable to be represented in proportional symbols
#' @param flagAbsSize Using the absolute value of the variable for graphic representation
#' @param flagCst Represent the location of the active samples only
#' @param useSel Use of the optional selection
#' @param asFactor Transform color variable into factor to use discrete palette
#' @param sizeRange Range of the size of the symbols (when 'nameSize' is not defined)
#' @param posX Rank of the coordinate used as the first coordinate
#' @param posY Rank of the coordinate used as the second coordinate
#' @param palette Name of the reference color map
#' @param naColor Color assigned to undefined samples
#' @param flagLegend Display the legend for point representation (color and size)
#' @param legendNameColor Name of the Legend for point representation by color
#' @param legendNameSize Name of the Legend for point representation by size
#' @param ... List of arguments passed to geom_point()
#' @return The description of the contents of the graphic layer
plot.symbol <- function(db, nameColor=NULL, nameSize=NULL, 
    flagAbsSize = FALSE, flagCst=FALSE, useSel=TRUE, asFactor=FALSE, 
    sizeRange=c(0.5, 3.), posX=0, posY=1,
    palette=NULL, naColor="transparent", flagLegend=FALSE, 
    legendNameColor=NA, legendNameSize=NA,  
    ...)
{ 
  p = list()

  # Get the name of the variable to be displayed by default
  nameSize = .defaultVariable(db, nameSize, nameColor)

  # Creating the necessary data frame
  df = .readPointCoor(db, useSel, posX, posY)
  aes_plt = aes(x=x, y=y)
  
  # Color of symbol
  if (! is.null(nameColor)) {
    colval  = db$getColumn(nameColor, TRUE)
    if (asFactor) colval = factor(colval)

    df["colour"] = colval
    aes_plt$colour = substitute(colval)
  }
  
  # Size of symbol
  if (! is.null(nameSize) && ! flagCst) {
    sizval  = db$getColumn(nameSize, TRUE)
    if (flagAbsSize) sizval = abs(sizval)

    df["size"] = sizval
    aes_plt$size = substitute(sizval)
  }

  layer <- geom_point(data = df, aes_plt, 
              na.rm=TRUE, show.legend = .showLegend(flagLegend), ...)
  p <- append(p, list(layer))

  # Define the Legend
  if (! is.null(nameColor)) {
    p <- append(p, .defineColour(palette, naColor = naColor,
      flagDiscrete = asFactor, title = legendNameColor))
    p <- .appendNewScale(p, "colour")
    p <- append(p, plot.decoration(title = nameColor, flagDefaultTitle = TRUE))
  }
  if (! is.null(nameSize) && ! flagCst) {
    p <- append(p, .defineSize(sizval, sizeRange, title = legendNameSize))
    p <- .appendNewScale(p, "size")
    p <- append(p, plot.decoration(title = nameSize, flagDefaultTitle = TRUE))
  }

  p
}

#' Plotting a point data base where samples are displayed with a label
#' @param db Data Base containing the information to be displayed
#' @param name Name of the variable to be represented
#' @param digit Number of decimal digits
#' @param useSel Use of the optional selection
#' @param posX Rank of the coordinate used as the first coordinate
#' @param posY Rank of the coordinate used as the second coordinate
#' @param ... List of arguments passed to geom_text_repel()
#' @return The description of the contents of the graphic layer
plot.literal <- function(db, name=NULL, digit=2, useSel=TRUE, posX=0, posY=1, 
    ...)
{
  p = list()

  if (!require(ggrepel, quietly=TRUE))
    stop("Package ggrepel is mandatory to use this function!")

  # Get the name of the variable to be displayed by default
  name = .defaultVariable(db, name)

  # Creating the necessary data frame
  df = .readPointCoor(db, useSel, posX, posY)
  aes_plt = aes(x=x, y=y)
  
  # Label of symbols
  labval  = round(db$getColumn(name,TRUE),digit)
  df["labval"] = as.character(labval)
  aes_plt$label = substitute(labval)
  
  layer <- geom_text_repel(data = df, aes_plt, na.rm=TRUE, ...)
  p <- append(p, list(layer))
  
  layer
}

#' Represent the contents of a variable defined on a grid as an Image
#' @param dbgrid Grid data base from gstlearn
#' @param name Name of the variable to be represented
#' @param useSel Use of an optional selection
#' @param posX rank of the coordinate which will serve as the first coordinate
#' @param posY rank of the coordinate which will serve as the second coordinate
#' @param corner A vector (same space dimension as 'dbgrid') which defines a pixel belonging to the extracted section
#' @param palette Name of the reference color map
#' @param naColor Color assigned to undefined samples
#' @param limits Bounds applied to the variable to be represented
#' @param legendName Name of the Legend for representation as an image
#' @param flagLegend Display the legend for grid representation as an image
#' @param ... Arguments passed to geom_raster() or geom_polygon()
#' @return The description of the contents of the figure
plot.raster <- function(dbgrid, name = NULL, useSel = TRUE, posX=0, posY=1, corner=NA, 
    palette=NULL, naColor="transparent", limits=NULL, legendName=NA, flagLegend=FALSE, ...)
{
  if (! .isGrid(dbgrid)) stop()

  if (!require(ggnewscale, quietly = TRUE)) 
    stop("Package ggnewscale is mandatory to use this function!")

  p = list()

  # Get the name of the variable to be displayed
  name = .defaultVariable(dbgrid, name)

  # Reading the Grid information
  df = .readGridCoor(dbgrid, name, useSel, posX, posY, corner)
  
  # Define the contents
  if (dbgrid$getAngles()[1] == 0 && ! dbgrid$hasSingleBlock())
  {
    layer <- geom_raster(
      data = df, mapping = aes(x = x, y = y, fill = data),
      show.legend = .showLegend(flagLegend), ...
    )
  }
  else
  {
    ids = seq(1, dbgrid$getNTotal())
    coords = dbgrid$getAllCellsEdges()
    positions = data.frame(id = rep(ids, each = 4), x = coords[[1]], y = coords[[2]])
    values = data.frame(id = ids, value = df$data)
    df <- merge(values, positions, by = c("id"))
    layer <- geom_polygon(
      data = df, mapping = aes(x = x, y = y, fill = value, group = id),
      show.legend = .showLegend(flagLegend), ...
    )
  }
  p <- append(p, list(layer))
  
  # Define the color Scale
  p <- append(p, .defineFill(palette, naColor = naColor, limits = limits, title = legendName))
  p <- .appendNewScale(p, "fill")
  
  p
}

#' Represent the contents of a variable defined on a grid with isovalues 
#' @param dbgrid Grid data base from gstlearn
#' @param name Name of the variable to be represented
#' @param useSel Use of an optional selection
#' @param posX rank of the coordinate which will serve as the first coordinate
#' @param posY rank of the coordinate which will serve as the second coordinate
#' @param corner A vector (same space dimension as 'dbgrid') which defines a pixel belonging to the extracted section
#' @param palette Name of the reference color map
#' @param naColor Color assigned to undefined samples
#' @param legendName Name of the Legend for representation as an image
#' @param flagLegend Display the legend for grid representation as an image
#' @param ... Arguments passed to geom_contour()
#' @return The description of the contents of the figure
plot.contour <- function(dbgrid, name=NULL, useSel = TRUE, posX=0, posY=1, corner=NA, 
    palette=NULL, naColor="transparent", legendName=NA, flagLegend=FALSE, ...)
{
  if (!.isGrid(dbgrid, TRUE)) stop()

  p = list()

  # Get the name of the variable to be displayed
  name = .defaultVariable(dbgrid, name)
    
  # Reading the Grid information
  df = .readGridCoor(dbgrid, name, useSel, posX, posY, corner)
  
  # Define the contents
  layer <- geom_contour(
    data = df, mapping = aes(x = x, y = y, z = data, colour=after_stat(level)), na.rm = TRUE,
    show.legend = .showLegend(flagLegend), ...
  )
  p <- append(p, list(layer))
  
  # Define the color Scale (only displayed if using various colors for)
  p <- append(p, .defineColour(palette, naColor = naColor, title = legendName))
  p <- .appendNewScale(p, "colour")

  p <- append(p, plot.decoration(title = name, flagDefaultTitle = TRUE))

  p
}

#' Representing a Polygon
#' @param poly A Polygon object from the gstlearn library
#' @param cols List of colors (optional)
#' @param flagTitle Plot the decoration attached to the figure
#' @param ... List of arguments passed to geom_polygon( )
#' @return The ggplot object
plot.polygon <- function(poly, cols=NA, flagTitle=FALSE, ...)
{
  dots = list(...)
  has_color = "color" %in% names(dots)
  
  p = list()
  npol = poly$getNPolyElem()
  if (missing(cols)) cols = .getColors()
  
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
#' @param ... List of arguments passed to geom_path() or geom_point()
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
#' @param flagLegend Show the legend when representing grid of occurrences (asPoint = FALSE)
#' @param legendName Name of the legend when representing grid of occurrences (asPoint = FALSE)
#' @param ... List of arguments passed to plot.XY() or plot.hist2d()
#' @return The ggplot object
plot.correlation <- function(db1, namex, namey, db2=NULL,
    asPoint=FALSE, 
    flagDiag=FALSE, diagColor = "red", diagLinetype = "solid", 
    flagRegr=FALSE, regrColor = "blue", regrLinetype = "solid", 
    flagBiss=FALSE, bissColor = "green", bissLinetype = "solid", 
    flagSameAxes=FALSE, flagLegend = FALSE, legendName = NA,
    ...)
{
  if (is.null(db2)) db2 = db1
  res = correlationPairs(db1, db2, namex, namey)
  x = db1$getValuesByNames(res[[1]], namex)
  y = db2$getValuesByNames(res[[2]], namey)
 
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
  if (flagLegend)
  {
    if (is.null(legendName)) legendName = "Count"
    p <- append(p, list(guides(fill = guide_colorbar(title=legendName, reverse=FALSE))))
  }
  else
  {
    p <- append(p, list(theme(legend.position='none')))
  }
  
  p 
}

#' Representing the H-scatter plot 
#' @param db A data base from gstlearn library
#' @param namex Name of the variable (within 'db1') which will be displayed along the horizontal axis
#' @param namey Name of the variable (within 'db2') which will be displayed along the vertical axis
#' @param varioparam A VarioParam structure describing the calculation criteria
#' @param ilag Rank of the lag to be used for calculations
#' @param idir Rank of the direction to be used for calculations
#' @param asPoint Represent samples pointwise if TRUE, otherwise as a grid painted with occurrences
#' @param flagDiag Represent the diagonal of the plot
#' @param diagColor Color of the diagonal
#' @param diagLinetype Line type of the diagonal
#' @param flagBiss Represent the first bisector (Y=X) 
#' @param bissColor Color of the first bisector (Y=X)
#' @param bissLinetype Line type of the first bisector (Y=X)
#' @param flagSameAxes Define the same bounds for horizontal and vertical axes
#' @param flagLegend Show the legend when representing grid of occurrences (asPoint = FALSE)
#' @param legendName Name of the legend when representing grid of occurrences (asPoint = FALSE)
#' @param ... List of arguments passed to plot.XY() or plot.hist2d()
#' @return The ggplot object
plot.hscatter <- function(db, namex, namey, varioparam, ilag=0, idir=0, asPoint=FALSE, 
    flagDiag=FALSE, diagColor = "red", diagLinetype = "solid", 
    flagBiss=FALSE, bissColor = "green", bissLinetype = "solid", 
    flagSameAxes=FALSE, flagLegend = FALSE, legendName = NA,
    ...)
{
  res = hscatterPairs(db, namex, namey, varioparam, ilag, idir)
  x = db$getValuesByNames(res[[1]], namex)
  y = db$getValuesByNames(res[[2]], namey)
 
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
  
  p = append(p, plot.decoration(xlab=namex, ylab=namey))
  
  # Set the Legend
  if (flagLegend)
  {
    if (is.null(legendName)) legendName = "Count"
    p <- append(p, list(guides(fill = guide_colorbar(title=legendName, reverse=FALSE))))
  }
  else
  {
    p <- append(p, list(theme(legend.position='none')))
  }
  
  p 
}

#' Representing a lithotype rule
#' @param rule A Rule object from gstlearn library
#' @param proportions The vector of facies proportions. When defined it is used to dimension the facies rectangles
#' @param maxG Maximum gaussian value (in absolute value)
#' @param cols List of colors (optional)
#' @param flagLegend Display the legend 
#' @param legendName Name of the Legend
#' @param ... List of arguments passed to geom_rect()
#' @return The ggplot object
plot.rule <- function(rule, proportions=NULL, maxG = 3., cols=NA, 
	flagLegend=FALSE, legendName="Facies", ...)
{
  p = list()
  nrect = rule$getNFacies()
  if (! is.null(proportions)) 
    rule$setProportions(proportions)
  else
    rule$setProportions()
  if (missing(cols)) cols = .getColors()
  
  df = data.frame(xmin=rep(0,nrect),xmax=rep(0,nrect),
      ymin=rep(0,nrect),ymax=rep(0,nrect), colors=rep(0,nrect))
  for (ifac in 1:nrect)
  {
    rect = rule$getThresh(ifac) # This function is 1-based
    df$xmin[ifac] = max(rect[1], -maxG)
    df$xmax[ifac] = min(rect[2], +maxG)
    df$ymin[ifac] = max(rect[3], -maxG)
    df$ymax[ifac] = min(rect[4], +maxG)
    df$colors[ifac] = ifac
  }
  
  p = append(p, geom_rect(data = df, mapping=aes(xmin = xmin, xmax = xmax, 
             ymin = ymin, ymax = ymax, fill = as.factor(colors)), na.rm=TRUE, ...))
     
  p <- append(p, .defineFill(cols, flagDiscrete = TRUE, ...))
  p <- .appendNewScale(p, "fill")  
  
  # Set the legend
  if (flagLegend)
  {
    p <- append(p, list(labs(color = legendName)))
  }
  else
  {
    p <- append(p, list(guides(color = "none")))
  }
  
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
plot.mesh <- function(mesh, flagFace=FALSE, flagApex=FALSE, rankMeshMax=-1, ...)
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
    if (flagCell && grid$isGrid())
    {
      edges = grid$getCellEdges(node)
      p = append(p, plot.XY(edges[[1]], edges[[2]], ...))
    }  
    
    # Represent the Neighborhood Ellipsoid
    if (neigh$getType()$getValue() == ENeigh_MOVING()$getValue())
    {
        edges = neigh$getEllipsoid(target)
        p = append(p, plot.XY(edges[[1]], edges[[2]], ...))
    
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
    }
    
    p
}

#' Represent the non-stationary parameters of a Model on a Grid
#' @param cova A Covariance object
#' @param dbgrid A Grid data base used for representation
#' @param useSel Use of an optional selection (masking off samples)
#' @param color Color used for graphic representation
#' @param flagOrtho Defines the long_axis of the anisotropy with respect to angle
#' @param scale Size given to the arraws
#' @param ... Arguments passed to geom_segment()
#' @return The ggplot object to geom_segment()
plot.covaOnGrid <- function(cova, dbgrid, useSel=TRUE, color='black', 
    flagOrtho=TRUE, scale=40, ...)
{
  if (!require(ggplot2, quietly = TRUE))
    stop("Package ggplot2 is mandatory to use this function!")
  
  # Extracting coordinates
  tab = dbgrid$getAllCoordinates(useSel)
  # Process the non-stationarity
  tabR1 = cova$informCoords(tab,EConsItem_RANGE(),0)
  tabR2 = cova$informCoords(tab,EConsItem_RANGE(),1)
  tabA  = cova$informCoords(tab,EConsItem_ANGLE())
  if (flagOrtho) tabA = 90 + tabA
  tabA = tabA * pi / 180.
  
  tabdx = (tabR1 * cos(tabA) - tabR2 * sin(tabA)) * scale
  tabdy = (tabR1 * sin(tabA) + tabR2 * cos(tabA)) * scale
  data = data.frame(x = tab[1,], y = tab[2,], dx=tabdx, dy=tabdy)
  #ax.quiver(tabx, taby, tabR2, tabR2, angles=tabA, color=color, **kwargs)
  
  p = ggplot(data = data, aes(x = x, y = y)) + 
      geom_point(size = 1) + 
      geom_segment(aes(xend = x + dx, yend = y + dy),
                   arrow = arrow(length = unit(0.1, "cm")), ...)
  
p
}

# The following Method definitions has been moved in rgstlearn.i (to prevent roxygen from crashing
# One example is provided next...
#setMethod("plot", signature(x="_p_AMesh"), function(x,y=missing,...)   plot.mesh(x,...))
