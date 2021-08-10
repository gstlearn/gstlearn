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

#include "AEnum.hpp"

/* Previously
typedef enum
{
  NEIGH_UNIQUE = 0,        //!< Unique Neighborhood
  NEIGH_BENCH = 1,         //!< Bench Neighborhood
  NEIGH_MOVING = 2,        //!< Moving Neighborhood
  NEIGH_IMAGE = 3,         //!< Image Neighborhood
} ENUM_NEIGHS;
*/

class ENeigh : public AEnum
{
  ENUM_DECLARE(ENeigh)

public:
  static ENeigh UNIQUE;        //!< Unique Neighborhood
  static ENeigh BENCH;         //!< Bench Neighborhood
  static ENeigh MOVING;        //!< Moving Neighborhood
  static ENeigh IMAGE;         //!< Image Neighborhood
};

