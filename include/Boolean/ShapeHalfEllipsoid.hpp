/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EShape.hpp"

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

  /// Interface for ICloneable
  IMPLEMENT_CLONING(ShapeHalfEllipsoid)

  EShape getType() const override { return EShape::fromKey("HALFELLIPSOID"); }
  int  getNParams() const override { return 4; }
  bool getFlagCutZ() const override { return true; }
  BooleanObject* generateObject(int ndim = 3) override;
  bool belongObject(const VectorDouble& coor, const BooleanObject* object) const override;
};
