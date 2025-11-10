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
#pragma once

#include "Transform/ATransform.hpp"
#include "Transform/ATransformWithAutoDiff.hpp"

#include "gstlearn_export.hpp"

namespace gstlrn
{

class GSTLEARN_EXPORT YeoJohnson:
#ifndef SWIG
  public ATransformWithAutoDiff<YeoJohnson>
#else
  public ATransform
#endif
{
public:
  YeoJohnson(double lambda);
  YeoJohnson(const YeoJohnson& r)            = delete;
  YeoJohnson& operator=(const YeoJohnson& r) = delete;
  virtual ~YeoJohnson()                      = default;

  double inverseTransform(double x) const override;
  double getLambdaValue() const { return _lambda.getValue(); }
  void setLambdaValue(double lambda) { _lambda.setValue(lambda); }
  void initParams() override;
  void appendParams(ListParams& listParams) override;

  String getName() const override { return "Yeo-Johnson"; }

  template<typename T>
  T evalImpl(T h) const
  {
    if (h >= 0)
    {
      if (_lambda.getValue() == 0)
      {
        return exp(h) - 1;
      }
      return pow( (_lambda.getValue() * h) + 1, 1 / _lambda.getValue()) - 1;
    }

    if (_lambda.getValue() == 2)
    {
      return 1 - exp(-h);
    }

    return 1 - pow(1 - ( (2 - _lambda.getValue()) * h ), 1 / (2 - _lambda.getValue()));
  }

private:
  ParamInfo _lambda;
};

} // namespace gstlrn