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

class GSTLEARN_EXPORT ShapeHalfSinusoid: public AShape
{
public:
  ShapeHalfSinusoid(double proportion = 1.,
                    double period = 10.,
                    double amplitude = 1.,
                    double thickness = 1.,
                    double xext = 1.,
                    double zext = 1.,
                    double theta = 0.);
  ShapeHalfSinusoid(const ShapeHalfSinusoid &r);
  ShapeHalfSinusoid& operator=(const ShapeHalfSinusoid &r);
  virtual ~ShapeHalfSinusoid();

  /// Interface for ICloneable
  IMPLEMENT_CLONING(ShapeHalfSinusoid)

  EShape getType() const override { return EShape::fromKey("HALFSINUSOID"); }
  int  getNParams() const override { return 6; }
  bool getFlagCutZ() const override { return true; }
  BooleanObject* generateObject(int ndim = 3) override;
  bool belongObject(const VectorDouble& coor, const BooleanObject* object) const override;
};
