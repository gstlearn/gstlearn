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
# This is a set of functions which enable performing 3D plots easily (from X. Freulon)
#

#' Initialize the visualization of the sphere using the DbGrid db containing
#' the apices of the triangles building the mesh.
#' The coordinates are:
#' x1 -> longitude in radiants (0 <= x1 <= 2*pi)
#' x2 -> colatitude in radiants(0 <= x1 <= pi)
#' @param val vector of values to be displayed
#' @value return a list with the function display
IniView_S2_in_3D <- function(grd) {
  stopifnot(grd$isGrid())
  stopifnot(grd$getNDim() == 2)
  
  # Getting the triangles and points
  mesh = MeshETurbo(grd)
  m1 = mesh$getAllApices()$toTL()
  m1_xyz = matrix(NaN, nrow = dim(m1)[1], ncol = 3)
  m1_xyz[,1] = sin(m1[,2]) * cos(m1[,1])
  m1_xyz[,2] = sin(m1[,2]) * sin(m1[,1])
  m1_xyz[,3] = cos(m1[,2])
  m2 = matrix(unlist(mesh$getAllMeshes()$getMatrix()), ncol = 3, byrow = TRUE)
  
  # Building the rgl mesh
  m3d = mesh3d(x = m1_xyz, triangles = 1 + t(m2))
  
  fn_display <- function(val, ncol = 50, palette = "heat") {
    display_on_S2(val = val, mesh = m3d, ncol = ncol, palette = palette)
  }
  list(display = fn_display)
}

#' mapping a real variable on the sphere
#' @param val vector of values to be displayed
#' @param mesh rgl on which the variable val is displayed
#' @param ncol number of colors used to display the variable
#' @param palette a string defining the palette
#' @value None but  the function rglwidget() should called outside
#' The number of values in val should be equal to the number of apices of mesh
display_on_S2 <- function(val, mesh, ncol = 50, palette = "heat") {
  stopifnot(class(mesh)[1] == "mesh3d")
  np = dim(mesh$vb)[2]
  stopifnot(length(val) == np)
  
  # Defining the colors
  zlim <- range(val)
  stopifnot(palette %in% c("heat", "rainbow", "terrain"))
  if (palette == "heat") {
    colorlut <- heat.colors(ncol) # height color lookup table
  } else if (palette == "rainbow") {
    colorlut <- rainbow(ncol) # height color lookup table
  } else if(palette == "terrain") {
    colorlut <- terrain.colors(ncol) # height color lookup table
  }
  
  # Assign colors to heights for each point
  col <- colorlut[ floor((ncol-1)*(val - zlim[1])/(zlim[2]-zlim[1])) + 1 ] 
  clear3d()
  open3d()
  shade3d(x = mesh, color = col)
  NULL  
}
