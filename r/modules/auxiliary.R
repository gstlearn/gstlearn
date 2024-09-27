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
# This is a set of functions which makes the life of gstlearn users easier
# They have been written by several authors and stored here
# Authors: Xavier FREULON
#

#' Plotting a series of experimental variograms
#' (to compare several simulations for example)
#'
#' @param var_list List containing the different experimental variograms
#' @param title    Title to be printed
#' @param idir     Rank of the direction to be plotted (0: all)
#'
compare_variograms <- function(var_list, title, idir = 0) {
  nsim  = length(var_list)
  ndir  = var_list[[1]]$getDirectionNumber()
  nlags = var_list[[1]]$getLagNumber(0)
  if (idir == 0) {
	ldir = 1:ndir
  } else {
        ldir = idir
  }

  # Loop over the directions
  for (idir in ldir) {
    # Get the lags
    hh = var_list[[1]]$getHhVec(idir = idir-1, ivar = 0, jvar = 0)
    # Get the Variograms for all variograms
    gg_sim = matrix(NaN, nrow = nlags-1, ncol = nsim)
    for (i in 1:nsim) {
      gg_sim[,i] = var_list[[i]]$getGgVec(idir = idir-1, ivar = 0, jvar = 0)
    }
    
    v_mean = apply(X = gg_sim, MARGIN = 1, FUN = mean)
    v_sd   = apply(X = gg_sim, MARGIN = 1, FUN = sd)
    
    # initial plot
    sdir = sprintf("Direction = %1d", idir)
    plot(NULL, NULL, xlim = range(hh), ylim = 1.2*c(0, 1.0),
         xlab = "Distance", ylab = "Variogram", xaxs="i", yaxs="i",
	 main = paste0(title, "\n", sdir))
    abline(h = 0, lty = 3, col = "gray")
    abline(v = 0, lty = 3, col = "gray")
    
    # experimental variogram
    lines(hh, v_mean, col = "orange", lw = 2)
    lines(hh, v_mean + 2 * v_sd, col = "orange", lty = 2)
    lines(hh, v_mean - 2 * v_sd, col = "orange", lty = 2)
    for (s in 1:nsim) {
      lines(hh, gg_sim[,s] , col = "gray", lty = 3)
    }
    
    # legend
    legend("bottomright", 
           legend = c("model", "mean variogram", "+/- 2 x Std.",
	   "empirical variograms"),
           lty = c(1, 1, 2, 3), lw = c(2, 2, 1, 1),
	   col = c("skyblue", "orange", "orange", "gray")
    )
  } # loop over the directions
  NULL
}

#' -----------------------------------------------------------------------------
#' Computing the variogram on a gridded sphere
#' -----------------------------------------------------------------------------
#' @param grd the grid
#' @param nm_var the name of the variable
#' @param nm_theta the name of the colatitude (in radians)
#' @param nm_phi the name of the longitude (in radians)
#' @param lag  a vector defining the *direction* in the theta-phi coordinate system
#' @param nlag number of lag to be computed
#' @param ndisc number of bins for the discretisation of the [0, dmax] interval
#' @param dmax maximum distance (in radians)
vario_on_S2 <- function(grd, nm_var, lag = c(1,0), nlag = 20, nm_theta = "theta", nm_phi = "phi", 
                        ndisc = 10, dmax = pi) {
  gg = rep(0.0, ndisc)
  hh = rep(0.0, ndisc)
  nc = rep(0, ndisc)
  stopifnot(grd$isGrid() & (grd$getNDim() == 2))
  nx = grd$getNXs()
  theta = matrix(grd[nm_theta], nrow = nx[1], ncol = nx[2])
  phi   = matrix(grd[nm_phi], nrow = nx[1], ncol = nx[2])
  var   = matrix(grd[nm_var], nrow = nx[1], ncol = nx[2])

  for (il in 1:nlag) {
    if((nx[1] > lag[1]*il)&(nx[2] > lag[2]*il)) {
      i1 = 1:(nx[1]-lag[1]*il)
      j1 = (1+lag[1]*il):nx[1]
      i2 = 1:(nx[2]-lag[2]*il)
      j2 = (1+lag[2]*il):nx[2]
      inc   = as.numeric(1/2 * (var[i1, i2] - var[j1, j2])^2)
      alpha = as.numeric(acos(sin(theta[i1,i2])*sin(theta[j1,j2])*cos(phi[i1,i2]-phi[j1,j2])+
                                cos(theta[i1,i2])*cos(theta[j1,j2])))
      for (n in 1:ndisc) {
        sel   = (alpha >= (n-1) * dmax/ndisc)&(alpha < (n)*dmax/ndisc)
        if(sum(sel)) {
          nc[n] = nc[n] + sum(sel)
          gg[n] = gg[n] + sum(inc[sel])
          hh[n] = hh[n] + sum(alpha[sel])
        }
      } # loop over the [0,pi] discretisation
    }
  } # loop over the lag
  sel = (nc > 0)
  gg[sel]  <- gg[sel]/nc[sel]
  hh[sel]  <- hh[sel]/nc[sel]

  list(gg = gg[sel], hh = hh[sel], nc = nc[sel], lag = lag, nlag = nlag, dmax = dmax, ndisc = ndisc)
}