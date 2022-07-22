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

/// TODO : use covariant return type to prevent dynamic_cast: 
// https://alfps.wordpress.com/2010/06/12/cppx-3-ways-to-mix-in-a-generic-cloning-implementation/
class GSTLEARN_EXPORT ICloneable
{
public:
  ICloneable() {};
  virtual ~ICloneable() {};

  virtual ICloneable* clone() const = 0;
};
