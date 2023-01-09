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

#include "Boolean/AShape.hpp"
#include "Boolean/ShapeParameter.hpp"

class BooleanObject;

class GSTLEARN_EXPORT ShapeParallelepiped: public AShape
{
public:
  ShapeParallelepiped(double proportion = 1.,
                      double xext = 1.,
                      double yext = 1.,
                      double zext = 1.,
                      double theta = 0.);
  ShapeParallelepiped(const ShapeParallelepiped &r);
  ShapeParallelepiped& operator=(const ShapeParallelepiped &r);
  virtual ~ShapeParallelepiped();

  /// Interface for ICloneable
  IMPLEMENT_CLONING(ShapeParallelepiped)

  EShape getType() const override { return EShape::fromKey("PARALLELEPIPED"); }
  int  getNParams() const override { return 4; }
  bool getFlagCutZ() const override { return false; }
  BooleanObject* generateObject(int ndim = 3) override;
  bool belongObject(const VectorDouble& coor, const BooleanObject* object) const override;
};
