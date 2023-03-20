/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EShape.hpp"

#include "Boolean/AShape.hpp"
#include "Boolean/ShapeParameter.hpp"

class BooleanObject;

class GSTLEARN_EXPORT ShapeEllipsoid: public AShape
{
public:
  ShapeEllipsoid(double proportion = 1.,
                 double xext = 1.,
                 double yext = 1.,
                 double zext = 1.,
                 double theta = 0.);
  ShapeEllipsoid(const ShapeEllipsoid &r);
  ShapeEllipsoid& operator=(const ShapeEllipsoid &r);
  virtual ~ShapeEllipsoid();

  /// Interface for ICloneable
  IMPLEMENT_CLONING(ShapeEllipsoid)

  EShape getType() const override { return EShape::fromKey("ELLIPSOID"); }
  int  getNParams() const override { return 4; }
  bool getFlagCutZ() const override { return false; }
  BooleanObject* generateObject(int ndim = 3) override;
  bool belongObject(const VectorDouble& coor, const BooleanObject* object) const override;
};
