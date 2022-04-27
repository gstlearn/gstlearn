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
#include "Boolean/AToken.hpp"
#include "Boolean/TokenParameter.hpp"

class GSTLEARN_EXPORT TokenParallelepiped: public AToken
{
public:
  TokenParallelepiped();
  TokenParallelepiped(const TokenParallelepiped &r);
  TokenParallelepiped& operator=(const TokenParallelepiped &r);
  virtual ~TokenParallelepiped();

  ETShape getType() const override { return ETShape::PARALLELEPIPED; }
  int getNArgs() const override { return 4; }
  int getFlagSymZ() const override { return 0; }
  String getParamName(int i) const override { return _paramNames[i]; }
  const TokenParameter& getParam(int i) const override { return _params[i]; }

private:
  VectorString _paramNames;
  std::vector<TokenParameter> _params;
};
