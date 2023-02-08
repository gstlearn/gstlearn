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

#include "Enum/EMorpho.hpp"

#include "geoslib_define.h"
#include "Arrays/BImage.hpp"
#include "Db/DbGrid.hpp"

GSTLEARN_EXPORT int morpho_count(const BImage& imagin);
GSTLEARN_EXPORT void morpho_duplicate(const BImage &imagin, BImage &imagout);
GSTLEARN_EXPORT VectorDouble morpho_labelling(int option,
                                              int flag_size,
                                              const BImage& imagin,
                                              double ccvoid,
                                              bool verbose = false);
GSTLEARN_EXPORT VectorInt morpho_labelsize(int option,
                                           const BImage& imagin);
GSTLEARN_EXPORT void morpho_erosion(int option,
                                    const VectorInt &radius,
                                    const BImage& imagin,
                                    BImage& imagout,
                                    bool verbose = false);
GSTLEARN_EXPORT void morpho_dilation(int option,
                                     const VectorInt &radius,
                                     const BImage& imagin,
                                     BImage& imagout,
                                     bool verbose = false);
GSTLEARN_EXPORT void morpho_intersection(const BImage& image1,
                                         const BImage& image2,
                                         BImage& imagout,
                                         bool verbose = false);
GSTLEARN_EXPORT void morpho_union(const BImage& image1,
                                  const BImage& image2,
                                  BImage &imagout,
                                  bool verbose = false);
GSTLEARN_EXPORT void morpho_opening(int option,
                                    const VectorInt &radius,
                                    const BImage& imagin,
                                    BImage& imagout);
GSTLEARN_EXPORT void morpho_closing(int option,
                                    const VectorInt &radius,
                                    const BImage& imagin,
                                    BImage& imagout);
GSTLEARN_EXPORT void morpho_negation(const BImage& imagin,
                                     BImage& imagout,
                                     bool verbse = false);
GSTLEARN_EXPORT void morpho_double2image(const VectorInt &nx,
                                         const VectorDouble &tab,
                                         double vmin,
                                         double vmax,
                                         BImage& imagout,
                                         bool verbose = false);
GSTLEARN_EXPORT BImage morpho_double2image(const VectorInt &nx,
                                           const VectorDouble &tab,
                                           double vmin,
                                           double vmax,
                                           bool verbose = false);
GSTLEARN_EXPORT void morpho_image2double(const BImage& imagin,
                                         int mode,
                                         double grain,
                                         double pore,
                                         VectorDouble &tab,
                                         bool verbose = false);

GSTLEARN_EXPORT void morpho_distance(int option,
                                     const VectorInt &radius,
                                     bool flag_dist_erode,
                                     BImage& imagin,
                                     VectorDouble &dist,
                                     bool verbose = false);
GSTLEARN_EXPORT void morpho_angle2D(DbGrid *dbgrid,
                                    const VectorInt &radius,
                                    int iptr0);
GSTLEARN_EXPORT void morpho_gradients(DbGrid *dbgrid, int iptr);

GSTLEARN_EXPORT VectorInt gridcell_neigh(int ndim,
                                         int option,
                                         int radius,
                                         int flag_center,
                                         bool verbose = false);

GSTLEARN_EXPORT int db_morpho_calc(DbGrid *dbgrid,
                                   int iptr0,
                                   const EMorpho& oper,
                                   double vmin = 0.,
                                   double vmax = 1.5,
                                   int option = 0,
                                   const VectorInt& radius = VectorInt(),
                                   bool dist_erode = false,
                                   bool verbose = false);
