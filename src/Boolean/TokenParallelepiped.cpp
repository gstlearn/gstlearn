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
#include "../../include/Boolean/TokenParallelepiped.hpp"

TokenParallelepiped::TokenParallelepiped()
    : AToken(),
      _paramNames(),
      _params()
{
  _paramNames.resize(4);
  _paramNames[0] = "X-Extension";
  _paramNames[1] = "Y-Extension";
  _paramNames[2] = "Z-Extension";
  _paramNames[3] = "Orientation Angle";
  _params.resize(4);
}

TokenParallelepiped::TokenParallelepiped(const TokenParallelepiped &r)
    : AToken(r),
      _paramNames(r._paramNames),
      _params(r._params)
{
}

TokenParallelepiped& TokenParallelepiped::operator=(const TokenParallelepiped &r)
{
  if (this != &r)
  {
    AToken::operator =(r);
    _paramNames = r._paramNames;
    _params = r._params;
  }
  return *this;
}

TokenParallelepiped::~TokenParallelepiped()
{
}

