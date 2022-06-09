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
#include "Basic/Vector.hpp"

class GSTLEARN_EXPORT IProjMatrix
{
public:
  IProjMatrix() { }
  virtual ~IProjMatrix() { }
  virtual int point2mesh(const VectorDouble& in, VectorDouble& out) const = 0;
  virtual int mesh2point(const VectorDouble& in, VectorDouble& out) const = 0;
  virtual int getApexNumber() const = 0;
  virtual int getPointNumber() const = 0;
};
