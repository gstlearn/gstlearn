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
#include "Boolean/AToken.hpp"
#include "Boolean/TokenParameter.hpp"

class Object;

class GSTLEARN_EXPORT TokenParaboloid: public AToken
{
public:
  TokenParaboloid(double proportion = 1.,
                  double xext = 1.,
                  double yext = 1.,
                  double zext = 1.,
                  double theta = 0.);
  TokenParaboloid(const TokenParaboloid &r);
  TokenParaboloid& operator=(const TokenParaboloid &r);
  virtual ~TokenParaboloid();

  /// Interface for Iclonable
  virtual IClonable* clone() const override { return new TokenParaboloid(*this); };

  ETShape getType() const override { return ETShape::PARABOLOID; }
  int  getNParams() const override { return 4; }
  bool getFlagCutZ() const override { return false; }
  Object* generateObject(int ndim = 3) override;
  bool belongObject(const VectorDouble& coor, const Object* object) const override;
};
