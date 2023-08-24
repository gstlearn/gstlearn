/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EMorpho.hpp"

#include "geoslib_define.h"
#include "Arrays/BImage.hpp"
#include "Db/DbGrid.hpp"

struct GSTLEARN_EXPORT Spill_Res
{
  bool   success; // TRUE if algorithm has been successfully performed
  double h;       // Elevation of the spill point
  double th;      // Maximum reservoir thickness
  int    ix0;     // Location of the Spill point grid node along X
  int    iy0;     // Location of the Spill point grid node along Y
};

/**
 * \defgroup MORPHO Mathematical Morphology:
 *
 * Here are the implementation of the traditional Mathematical Morphology functions.
 *
 * These operations can be performed using Binary Image format (BImage) or colored arrays (double).
 * Some functions are applied to DbGrid (the input if a field of the DbGrid and results are written as
 * new fields in the same DbGrid); most functions are activated on BImage directly.
 *
 * For examples, see #!TUTORIALS!Tuto_Morpho
 *
 **/

/** @addtogroup MORPHO_0 Initial manipulations
 * \ingroup MORPHO
 *
 * @param  imagin     Pointer to the BImage containing the input image
 * @param  imagout    Pointer to the BImage which will receive the output image
 *                    (It must have been allocated beforehand: same dimension as 'imagin')
 *  @{
 */
GSTLEARN_EXPORT int morpho_count(const BImage& imagin);
GSTLEARN_EXPORT void morpho_duplicate(const BImage &imagin, BImage &imagout);
/**@}*/

/** @addtogroup MORPHO_1 Basic operations
 * \ingroup MORPHO
 *
 * @param  imagin     Pointer to the BImage containing the input image
 * @param  imagout    Pointer to the BImage which will receive the output image
 *                    (It must have been allocated beforehand: same dimension as 'imagin')
 * @param  option     Description of the structuring element:
 *                    1 for BLOCK and 0 for CROSS (see remarks for more information)
 * @param  radius     Vector giving the extensions of the structuring element along each direction
 *                    of the Space.
 * @param  verbose    Verbose flag
 *
 * @remarks In 2-D, if the central cell is denoted (i,j) and the radius is set to 1, the structuring element contains:
 *          - for CROSS: (i+1,j); (i-1,j); (i,j+1); (i,j-1)
 *          - for BLOCK: (i-1,j-1); (i-1,j); (i-1,j+1); (i,j-1); (i,j+1); (i+1,j-1); (i+1,j); (i+1,j+1)
 *  @{
 */
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
GSTLEARN_EXPORT void morpho_opening(int option,
                                    const VectorInt &radius,
                                    const BImage& imagin,
                                    BImage& imagout,
                                    bool verbose = false);
GSTLEARN_EXPORT void morpho_closing(int option,
                                    const VectorInt &radius,
                                    const BImage& imagin,
                                    BImage& imagout,
                                    bool verbose = false);
/**@}*/

/** @addtogroup MORPHO_2 Logical operations on images
 * \ingroup MORPHO
 *
 * @param  image1     Pointer to the BImage containing the first input image (for diadic operation)
 * @param  image2     Pointer to the BImage containing the second input image (for diadic opertion)
 * @param  imagin     Pointer to the BImage containing the one input image (for monadic operation)
 * @param  imagout    Pointer to the BImage which will receive the output image
 *                    (It must have been allocated beforehand: same dimension as 'imagin')
 * @param  verbose    Verbose flag
 *  @{
 */
GSTLEARN_EXPORT void morpho_intersection(const BImage& image1,
                                         const BImage& image2,
                                         BImage& imagout,
                                         bool verbose = false);
GSTLEARN_EXPORT void morpho_union(const BImage& image1,
                                  const BImage& image2,
                                  BImage &imagout,
                                  bool verbose = false);
GSTLEARN_EXPORT void morpho_negation(const BImage& imagin,
                                     BImage& imagout,
                                     bool verbose = false);
/**@}*/

/** @addtogroup MORPHO_3 Convert arrays (provided as a Vector of double values) into Bimage (binary image) and vice-versa
 * \ingroup MORPHO
 *
 * @param  nx         Vector giving the number of nodes for each space dimension
 * @param  tabin      Array of double values containing the input information
 * @param  tabout     Array of double values which will receive the output values
 *                    (It must have been allocated beforehand)
 * @param  imagin     Pointer to the BImage containing the one input image
 * @param  imagout    Pointer to the BImage which will receive the output image.
 *                    This is used when converting "to" a BImage (in place)
 *                    (It must have been allocated beforehand: same dimension as 'imagin')
 *
 * @param  vmin, vmax If the value of 'tab' lies within [vmin, vmax], BImage is set to 1; 0 otherwise
 *                    (used when converting from double to BImage)
 * @param  grain, pore Double values assigned to 1 or 0 of the input BImage
 *                     (used when converting from BImage to double)
 * @param  verbose    Verbose flag
 *  @{
 */
GSTLEARN_EXPORT void morpho_double2imageInPlace(const VectorInt &nx,
                                                const VectorDouble &tabin,
                                                double vmin,
                                                double vmax,
                                                BImage &imagout,
                                                bool verbose = false);
GSTLEARN_EXPORT BImage morpho_double2image(const VectorInt &nx,
                                           const VectorDouble &tabin,
                                           double vmin,
                                           double vmax,
                                           bool verbose = false);
GSTLEARN_EXPORT void morpho_image2double(const BImage& imagin,
                                         int mode,
                                         double grain,
                                         double pore,
                                         VectorDouble &tabout,
                                         bool verbose = false);
/**@}*/

/** @addtogroup MORPHO_4 Label the connected component of a BImage
 * \ingroup MORPHO
 *
 * @param  imagin     Pointer to the BImage containing the one input image
 * @param  option     Description of the structuring element:
 *                    1 for BLOCK and 0 for CROSS (see remarks for more information)
 * @param  flag_size  1 if output information is labeled in volume of the connected component
 *                    0 if output information contains the rank of the connected component
 * @param  ccvoid     Value assigned to pixels which do not belong to any
 *                    connected component
 * @param  verbose    Verbose flag
 *
 * @remarks In 2-D, if the central cell is denoted (i,j) and the radius is set to 1, the structuring element contains:
 *          - for CROSS: (i+1,j); (i-1,j); (i,j+1); (i,j-1)
 *          - for BLOCK: (i-1,j-1); (i-1,j); (i-1,j+1); (i,j-1); (i,j+1); (i+1,j-1); (i+1,j); (i+1,j+1)
 *
 *  @{
 */
GSTLEARN_EXPORT VectorDouble morpho_labelling(int option,
                                              int flag_size,
                                              const BImage& imagin,
                                              double ccvoid,
                                              bool verbose = false);
GSTLEARN_EXPORT VectorInt morpho_labelsize(int option,
                                           const BImage& imagin);
/**@}*/

/** @addtogroup MORPHO_5 Performs some miscellaneous operations on BImages
 * \ingroup MORPHO
 *
 * @param  imagin     Pointer to the BImage containing the one input image
 * @param  option     Description of the structuring element:
 *                    1 for BLOCK and 0 for CROSS (see remarks for more information)
 * @param  radius     Vector giving the extensions of the structuring element along each direction
 *                    of the Space.
 * @param flagDistErode Inflate the grain; 0 Reduce the grain
 * @param dist Contains the vector of returned distances
 * @param iptr0 UID where the angle calculation will be calculated
 * @param flag_center True to omit the center
 *
 * @param  verbose    Verbose flag
 *
 * @remarks In 2-D, if the central cell is denoted (i,j) and the radius is set to 1, the structuring element contains:
 *          - for CROSS: (i+1,j); (i-1,j); (i,j+1); (i,j-1)
 *          - for BLOCK: (i-1,j-1); (i-1,j); (i-1,j+1); (i,j-1); (i,j+1); (i+1,j-1); (i+1,j); (i+1,j+1)
 *
 *  @{
 */
GSTLEARN_EXPORT void morpho_distance(int option,
                                     const VectorInt &radius,
                                     bool flagDistErode,
                                     BImage& imagin,
                                     VectorDouble &dist,
                                     bool verbose = false);
GSTLEARN_EXPORT VectorInt gridcell_neigh(int ndim,
                                         int option,
                                         int radius,
                                         bool flag_center = true,
                                         bool verbose = false);
/**@}*/

GSTLEARN_EXPORT Spill_Res spillPoint(DbGrid *dbgrid,
                                     const String& name_depth,
                                     const String& name_data,
                                     int option = 0,
                                     bool flag_up = true,
                                     int verbose_step = 0,
                                     double hmax = TEST);
