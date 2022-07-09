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

class GSTLEARN_EXPORT ShapeHalfEllipsoid: public AShape
{
public:
  ShapeHalfEllipsoid(double proportion = 1.,
                     double xext = 1.,
                     double yext = 1.,
                     double zext = 1.,
                     double theta = 0.);
  ShapeHalfEllipsoid(const ShapeHalfEllipsoid &r);
  ShapeHalfEllipsoid& operator=(const ShapeHalfEllipsoid &r);
  virtual ~ShapeHalfEllipsoid();

  /// Interface for IClonable
  virtual IClonable* clone() const override { return new ShapeHalfEllipsoid(*this); };

  ETShape getType() const override { return ETShape::HALFELLIPSOID; }
  int  getNParams() const override { return 4; }
  bool getFlagCutZ() const override { return true; }
  BooleanObject* generateObject(int ndim = 3) override;
  bool belongObject(const VectorDouble& coor, const BooleanObject* object) const override;
};
