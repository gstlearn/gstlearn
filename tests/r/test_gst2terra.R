#' test_gst2terra.R
#' To test the consistency between the initial DbGrid and raster grid
#' two functions are tested: *gstlearn_2_terra* and *terra_to_gstlearn*
#' Author: X. FREULON
#' Creation: 2025/06/05

rm(list = ls())
suppressWarnings(suppressMessages(library(terra)))
suppressWarnings(suppressMessages(library(sf)))
suppressWarnings(suppressMessages(library(gstlearn)))

# creation of the initial 2d DbGrid
gst_ini = gstlearn::DbGrid_create(nx=c(7,11), dx=c(2,5), x0 = c(100,-2588))
gst_ini["Simu"] = runif(n = gst_ini$getNSample())
# creation of the initial terra raster
terra_ini = gstlearn::gstlearn_to_terra(gst_ini)
# creation of the final gstlearn grid
gst_final = gstlearn::terra_to_gstlearn(terra_ini)
# creation of the final terra raster
terra_final = gstlearn::gstlearn_to_terra(gst_final)

# comparison of the DbGrid paramaters
stopifnot(all(gst_ini$getX0s() == gst_final$getX0s()))
stopifnot(all(gst_ini$getNXs() == gst_final$getNXs()))
stopifnot(all(gst_ini$getDXs() == gst_final$getDXs()))

# comparison of the raster paramaters
stopifnot(all(st_bbox(terra_ini) == st_bbox(terra_final)))
stopifnot(all(dim(terra_ini) == dim(terra_final)))
stopifnot(all(res(terra_ini) == res(terra_final)))

# consistency between DbGrid and raster
stopifnot(all(gst_ini$getNXs() == dim(terra_final)[c(2,1)]))
stopifnot(all(gst_ini$getDXs() == res(terra_final)))

print("Test is successful!")