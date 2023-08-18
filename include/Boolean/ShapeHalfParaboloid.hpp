/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EShape.hpp"

#include "Boolean/AShape.hpp"
#include "Boolean/ShapeParameter.hpp"

class BooleanObject;

class GSTLEARN_EXPORT ShapeHalfParaboloid: public AShape
{
public:
  ShapeHalfParaboloid(double proportion = 1.,
                      double xext = 1.,
                      double yext = 1.,
                      double zext = 1.,
                      double theta = 0.);
  ShapeHalfParaboloid(const ShapeHalfParaboloid &r);
  ShapeHalfParaboloid& operator=(const ShapeHalfParaboloid &r);
  virtual ~ShapeHalfParaboloid();

  /// Interface for ICloneable
  IMPLEMENT_CLONING(ShapeHalfParaboloid)

  EShape getType() const override { return EShape::fromKey("HALFPARABOLOID"); }
  int  getNParams() const override { return 4; }
  bool getFlagCutZ() const override { return true; }
  BooleanObject* generateObject(int ndim = 3) override;
  bool belongObject(const VectorDouble& coor, const BooleanObject* object) const override;
};
