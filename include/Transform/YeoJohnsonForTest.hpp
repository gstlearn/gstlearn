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
#include "gstlearn_export.hpp"

#include <cmath>
// The aim of this class is to test the generic mechanisms of ATransform and to compare their result
// with the analytical version implemented in the YeoJohnson official class.
namespace gstlrn
{

class YeoJohnsonForTest: public ATransform
{
public:
  YeoJohnsonForTest(double lambda)
    : _lambda(lambda) {};
  YeoJohnsonForTest(const YeoJohnsonForTest& r)            = default;
  YeoJohnsonForTest& operator=(const YeoJohnsonForTest& r) = default;
  IMPLEMENT_CLONING(YeoJohnsonForTest)
  String getName() const override { return "YeoJohnsonForTest"; }
  virtual ~YeoJohnsonForTest() = default;

  bool hasParameters() const override { return false; }
  double getLambda() const { return _lambda; }
  void setLambda(double lambda) { _lambda = lambda; }

  double transform(double h) const override
  {
    if (h >= 0)
    {
      if (_lambda == 0)
      {
        return exp(h) - 1;
      }
      return pow((_lambda * h) + 1, 1 / _lambda) - 1;
    }

    if (_lambda == 2)
    {
      return 1 - exp(-h);
    }

    return 1 - pow(1 - ((2 - _lambda) * h), 1 / (2 - _lambda));
  }

private:
  double _lambda;
};

} // namespace gstlrn
