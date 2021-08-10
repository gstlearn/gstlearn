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

#include "Enum/ENeigh.hpp"

/* Previously
typedef enum
{
  NEIGH_UNIQUE = 0,        //!< Unique Neighborhood
  NEIGH_BENCH = 1,         //!< Bench Neighborhood
  NEIGH_MOVING = 2,        //!< Moving Neighborhood
  NEIGH_IMAGE = 3,         //!< Image Neighborhood
} ENUM_NEIGHS;
*/

ENUM_INIT(ENeigh)

ENUM_DEFINE(ENeigh, UNIQUE, 0, "Unique neighborhood")
ENUM_DEFINE(ENeigh, BENCH , 1, "Bench neighborhood")
ENUM_DEFINE(ENeigh, MOVING, 2, "Moving neighborhood")
ENUM_DEFINE(ENeigh, IMAGE , 3, "Image neighborhood")

ENUM_DEFAULT(ENeigh, UNIQUE)
