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
#include "Boolean/ETShape.hpp"
#include "Boolean/AShape.hpp"
#include "Boolean/ShapeParameter.hpp"

class BooleanObject;

class GSTLEARN_EXPORT ShapeParaboloid: public AShape
{
public:
  ShapeParaboloid(double proportion = 1.,
                  double xext = 1.,
                  double yext = 1.,
                  double zext = 1.,
                  double theta = 0.);
  ShapeParaboloid(const ShapeParaboloid &r);
  ShapeParaboloid& operator=(const ShapeParaboloid &r);
  virtual ~ShapeParaboloid();

  /// Interface for IClonable
  virtual IClonable* clone() const override { return new ShapeParaboloid(*this); };

  ETShape getType() const override { return ETShape::PARABOLOID; }
  int  getNParams() const override { return 4; }
  bool getFlagCutZ() const override { return false; }
  BooleanObject* generateObject(int ndim = 3) override;
  bool belongObject(const VectorDouble& coor, const BooleanObject* object) const override;
};
