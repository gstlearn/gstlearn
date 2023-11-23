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
# This is a set of functions for interface with between gstlean package
# and geomatic package. In particular if converts back-and-forth
# information from
# - sf
# - terra
# Author: Xavier FREULON

#' Convert a *sf* object into a *gstlearn* object. 
#' The POLYGON and MULTIPOLYGON are converted into a Polygons.
#' The crs information and attributes are not stored.
#' The POINT and MULTIPOINT are converted into a Db.
#' The crs information is not stored.
#'
#' @param x a sf object to be converted
#' @param quiet a Boolean to control messages
#'
#' @value returns the gstlearn object (Polygons or Db)
sf_to_gstlearn <- function(x, quiet = TRUE)
{
  if (!require(sf, quietly=TRUE))
    stop("Package 'sf' is mandatory to use this function!")

  val = NULL
  
  # Conversion of polygons
  geo = st_geometry(x)
  if(class(geo)[1] == "sfc_MULTIPOLYGON") {
    if (! quiet) print(paste("Converting polygons..."))
    xy  = st_coordinates(geo)
    idx = sort(unique(xy[,4]))
    val = Polygons()
    for (i in idx) {
      sel <- (xy[, 4] == i)
      polyElem <- PolyElem(x = xy[sel, 1], y = xy[sel, 2])
      err = val$addPolyElem(polyElem)
      }
    }

  # conversion of scattered points
  if((class(geo)[1] == "sfc_POINT")) {
    if (! quiet) print(paste("Converting scattered points..."))
    val = Db()
    xy = st_coordinates(x)
    val[colnames(xy)[1]] <- xy[,1]
    val[colnames(xy)[2]] <- xy[,2]
    err = val$setLocators(names = colnames(xy), locatorType = ELoc_X(),
          cleanSameLocator = TRUE)
    dat = st_drop_geometry(x)
    vn = colnames(dat)
    for (v in vn) {
      values <- unlist(st_drop_geometry(dat[v]))
      if(is.numeric(values)) {
        val[v] <- values
      } else {
        print(paste(">>> variable ", v, " is not converted (not numeric)"))  
      }
    }
  }
  val
}

#' Convert a *gstlearn* object into a *sf* object. 
#' The Polygons is converted into a MULTIPOLYGON. 
#' The Db is converted into a MULTIPOINT.
#' The crs information should be provided.
#'
#' @param x a gstlearn object to be converted
#' @param crs a string storing the reference to the coordinates system
#' (e.g. "EPSG:4326" for long/lat in WGS84)
#'
#' @value returns the sf object
gstlearn_to_sf <- function(x, crs = NA)
{
  if (!require(sf, quietly=TRUE))
    stop("Package 'sf' is mandatory to use this function!")

  val = NULL
  if (class(x)[1] == "_p_Polygons") {
      lp = list()
      for (i in 0:(x$getPolyElemNumber()-1)) {
        lp[[1+length(lp)]] <- list(matrix(c(x$getX(i), x$getY(i)),
	                           ncol = 2, byrow = FALSE))
      }
      val = st_sf(geometry = st_sfc(st_multipolygon(lp))) |>
          st_set_crs(crs)
  }
  
  if (class(x)[1] == "_p_Db") {
    df <- x[]
    val <- st_as_sf(df, coords = x$getNamesByLocator(ELoc_X())) |>
      st_set_crs(crs)
  }
  if (is.null(val)) {
   print(paste("gstlearn_to_sf: class = ", class(x)[1], " not yet implemented."))
 }
 val
}

#' Convert a *terra* raster into a *gstlearn* DbGrid. 
#' The crs information is not stored.
#'
#' @param x a terra raster to be converted
#'
#' @value returns the DbGrid object
terra_to_gstlearn <- function(x)
{
  if (!require(terra, quietly=TRUE))
    stop("Package 'terra' is mandatory to use this function!")

  stopifnot(class(x)[1] == "SpatRaster")
  dx = res(x)[c(2,1)]
  nx = dim(x)[c(2,1)]
  x0 = st_bbox(x)[1:2]
  grd = DbGrid_create(nx = nx, dx = dx, x0 = x0)
  for (i in 1:dim(x)[3]) {
    val_ini = x[[names(x)[i]]][]
    val = rep(NaN, length(val_ini))
    for (j in 1:nx[2]) {
      idx_ini = 1:nx[1] + nx[1]*(j - 1)
      idx     = 1:nx[1] + nx[1]*(nx[2] - j) 
      val[idx] <- val_ini[idx_ini]
    }
    grd[names(x)[i]] = val
  }
 grd
}

#' Convert a *gstlearn* grid into a *terra* raster. (2 dimensional rasters) 
#' The crs information should be provided.
#'
#' @param x a gstlearn DbGrid to be converted
#' @param crs a string storing the reference to the coordinates system
#' (e.g. "EPSG:4326" for long/lat in WGS84)
#'
#' @value returns the sf object
gstlearn_to_terra <- function(x, crs = NA)
{
  if (!require(terra, quietly=TRUE))
    stop("Package 'terra' is mandatory to use this function!")

  stopifnot(x$isGrid() & (x$getNDim() == 2))
  nx    <- x$getNXs()
  dx    <- x$getDXs()
  x_min <- x$getX0s()
  x_max <- x_min + dx * (nx - 1)
  r <- rast(nrows = nx[2] , ncols = nx[1], 
            xmin = x_min[1], xmax = x_max[1], ymin = x_min[2], ymax = x_max[2])
  crs(r) <- crs
  var_nm <- x$getAllNames()[-c(1,2,3)]
  for (v in var_nm) {
    val_ini = x[v]
    val = rep(NaN, length(val_ini))
    for (j in 1:nx[2]) {
      idx_ini = 1:nx[1] + nx[1]*(j - 1)
      idx     = 1:nx[1] + nx[1]*(nx[2] - j) 
      val[idx] <- val_ini[idx_ini]
      }
    r[[v]] <- val
    }
  r
}
