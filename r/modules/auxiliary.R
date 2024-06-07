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
compare_variograms <- function(var_list, title, dir = 0) {
  nsim = length(var_list)
  ndir = var_list[[0]]$getDirectionNumber()
  nlags = var_list[[0]]$getLagNumber(idir-1)
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
