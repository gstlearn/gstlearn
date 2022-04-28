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

class GSTLEARN_EXPORT TokenHalfEllipsoid: public AToken
{
public:
  TokenHalfEllipsoid();
  TokenHalfEllipsoid(const TokenHalfEllipsoid &r);
  TokenHalfEllipsoid& operator=(const TokenHalfEllipsoid &r);
  virtual ~TokenHalfEllipsoid();

  ETShape getType() const override { return ETShape::HALFELLIPSOID; }
  int  getNParams() const override { return 4; }
  bool getFlagCutZ() const override { return true; }
  Object* generateObject(int ndim = 3) override;
  bool belongObject(const VectorDouble& coor, const Object* object) const override;
};
