/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/

#include "Transform/YeoJohnson.hpp"
#include "Basic/ParamInfo.hpp"
#include "Transform/ATransform.hpp"
#include "geoslib_define.h"

namespace gstlrn
{

YeoJohnson::YeoJohnson(double lambda)
  : ATransformWithAutoDiff<YeoJohnson>()
  , _lambda(ParamInfo("Yeo-Johnson Lambda", lambda,{-1., 2.}))
{
}

void YeoJohnson::initParams()
{
  _lambda.setValue(0.0);
}

double YeoJohnson::inverseTransform(double x) const
{
  if (x >= 0)
  {
    if (_lambda.getValue() == 0)
    {
      return std::log(x + 1);
    }
    return (std::pow(x + 1, _lambda.getValue()) - 1) / _lambda.getValue();
  }

  if (_lambda.getValue() == 2)
  {
    return -std::log(-x + 1);
  }

  return -(std::pow(-x + 1, 2 - _lambda.getValue()) - 1) / (2 - _lambda.getValue());
}

void YeoJohnson::appendParams(ListParams& listParams)
{
  listParams.addParam(_lambda);
}

} // namespace gstlrn
