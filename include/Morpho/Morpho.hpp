/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

GSTLEARN_EXPORT int morpho_image_size(const VectorInt &nx);
GSTLEARN_EXPORT int morpho_count(const VectorInt &nx,
                                 const VectorUChar &imagin);
GSTLEARN_EXPORT VectorUChar morpho_image_manage(const VectorInt &nx);
GSTLEARN_EXPORT void morpho_duplicate(const VectorInt &nx,
                                      const VectorUChar &imagin,
                                      VectorUChar &imagout);
GSTLEARN_EXPORT int morpho_labelling(const VectorInt &nx,
                                     int option,
                                     int flag_size,
                                     const VectorUChar &imagin,
                                     double ccvoid,
                                     VectorDouble &compnum,
                                     bool verbose = false);
GSTLEARN_EXPORT VectorDouble morpho_labelling(const VectorInt &nx,
                                              int option,
                                              int flag_size,
                                              const VectorUChar &imagin,
                                              double ccvoid,
                                              bool verbose = false);
GSTLEARN_EXPORT VectorInt morpho_labelsize(const VectorInt &nx,
                                           int option,
                                           const VectorUChar &imagin);
GSTLEARN_EXPORT void morpho_erosion(const VectorInt &nx,
                                    int option,
                                    const VectorInt &radius,
                                    const VectorUChar &imagin,
                                    VectorUChar &imagout,
                                    bool verbose = false);
GSTLEARN_EXPORT void morpho_dilation(const VectorInt &nx,
                                     int option,
                                     const VectorInt &radius,
                                     const VectorUChar &imagin,
                                     VectorUChar &imagout,
                                     bool verbose = false);
GSTLEARN_EXPORT void morpho_intersection(const VectorInt &nx,
                                         const VectorUChar &image1,
                                         const VectorUChar &image2,
                                         VectorUChar &imagout,
                                         bool verbose = false);
GSTLEARN_EXPORT void morpho_union(const VectorInt &nx,
                                  const VectorUChar &image1,
                                  const VectorUChar &image2,
                                  VectorUChar &imagout,
                                  bool verbose = false);
GSTLEARN_EXPORT void morpho_opening(const VectorInt &nx,
                                    int option,
                                    const VectorInt &radius,
                                    const VectorUChar &imagin,
                                    VectorUChar &imagout);
GSTLEARN_EXPORT void morpho_closing(const VectorInt &nx,
                                    int option,
                                    const VectorInt &radius,
                                    const VectorUChar &imagin,
                                    VectorUChar &imagout);
GSTLEARN_EXPORT void morpho_negation(const VectorInt &nx,
                                     const VectorUChar &imagin,
                                     VectorUChar &imagout,
                                     bool verbse = false);
GSTLEARN_EXPORT void morpho_double2image(const VectorInt &nx,
                                         const VectorDouble &tab,
                                         double vmin,
                                         double vmax,
                                         VectorUChar &imagout,
                                         bool verbose = false);
GSTLEARN_EXPORT VectorUChar morpho_double2image(const VectorInt &nx,
                                                const VectorDouble &tab,
                                                double vmin,
                                                double vmax,
                                                bool verbose = false);
GSTLEARN_EXPORT void morpho_image2double(const VectorInt &nx,
                                         const VectorUChar &imagin,
                                         int mode,
                                         double grain,
                                         double pore,
                                         VectorDouble &tab,
                                         bool verbose = false);
GSTLEARN_EXPORT void morpho_distance(const VectorInt &nx,
                                     int option,
                                     const VectorInt &radius,
                                     int flag_erode,
                                     VectorUChar &imagin,
                                     VectorDouble &dist);
GSTLEARN_EXPORT void morpho_angle(const VectorInt &nx,
                                  int radius,
                                  double *tab,
                                  double *tabout);
GSTLEARN_EXPORT void bitmap_print(const VectorInt &nx,
                                  const VectorUChar &imagin);
GSTLEARN_EXPORT int bitmap_size(const VectorInt &nx);
GSTLEARN_EXPORT int bitmap_get_value(const VectorInt &nx,
                                     const VectorUChar &imagin,
                                     int ix,
                                     int iy,
                                     int iz);
GSTLEARN_EXPORT void bitmap_set_value(const VectorInt &nx,
                                      VectorUChar &imagout,
                                      int ix,
                                      int iy,
                                      int iz,
                                      int bitval);
GSTLEARN_EXPORT VectorInt gridcell_neigh(int ndim,
                                         int option,
                                         int radius,
                                         int flag_center,
                                         int verbose,
                                         int *nvois);
