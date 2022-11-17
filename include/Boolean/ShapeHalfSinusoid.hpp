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

#include "Enum/ETShape.hpp"

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

  ETShape getType() const override { return ETShape::HALFSINUSOID; }
  int  getNParams() const override { return 6; }
  bool getFlagCutZ() const override { return true; }
  BooleanObject* generateObject(int ndim = 3) override;
  bool belongObject(const VectorDouble& coor, const BooleanObject* object) const override;
};
