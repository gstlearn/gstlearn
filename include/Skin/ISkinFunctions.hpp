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

class GSTLEARN_EXPORT ISkinFunctions
{
public:
  ISkinFunctions() {};
  virtual ~ISkinFunctions() {};

  virtual int isAlreadyFilled(int /*ipos*/) const = 0;
  virtual int isToBeFilled(int /*ipos*/) const = 0;
  virtual double getWeight(int /*ipos*/, int /*idir*/) const { return 1; }
};
