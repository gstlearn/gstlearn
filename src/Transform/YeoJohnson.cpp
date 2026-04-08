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
    , _lambda(ParamInfo("Yeo-Johnson Lambda", lambda, {-1., 2.}))
  {
  }

  YeoJohnson::YeoJohnson(const YeoJohnson& r)
    : ATransformWithAutoDiff<YeoJohnson>(r)
    , _lambda(r._lambda)
    , _K(r._K)
  {
  }

  YeoJohnson& YeoJohnson::operator=(const YeoJohnson& r)
  {
    if (this != &r)
    {
      ATransformWithAutoDiff<YeoJohnson>::operator=(r);
      _lambda = r._lambda;
      _K = r._K;
    }
    return *this;
  }

  YeoJohnson* YeoJohnson::create(double lambda)
  {
    return new YeoJohnson(lambda);
  }

  void YeoJohnson::initParams(double min, double max)
  {
    DECLARE_UNUSED(min);
    if (max != INF) _K = 1.2 * max;
  }

  void YeoJohnson::_printParams(
    std::stringstream& sstr,
    const AStringFormat* strfmt) const
  {
    DECLARE_UNUSED(strfmt);
    sstr << "Parameters:\n";
    sstr << "  Lambda: " << _lambda.getValue() << "\n";
    sstr << "  Saturation: " << _K << "\n";
  }

  VectorDouble YeoJohnson::getParams() const
  {
    return {_lambda.getValue()};
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

    return -(std::pow(-x + 1, 2 - _lambda.getValue()) - 1)
         / (2 - _lambda.getValue());
  }

  void YeoJohnson::appendParams(ListParams& listParams)
  {
    listParams.addParam(_lambda);
  }

} // namespace gstlrn
