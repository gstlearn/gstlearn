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

#include "Basic/ParamInfo.hpp"
#include "Transform/ATransform.hpp"
#include "Transform/ATransformWithAutoDiff.hpp"

#include "gstlearn_export.hpp"

// This class implements the Tuckey GH transformation
// See: Xu, G., & Genton, M. G. (2017) Tuckey g-and-h Random fields.
// Journal of the American Statistical Association, 112 (519), 1236-1249.

namespace gstlrn
{

class GSTLEARN_EXPORT TuckeyGH:
#ifndef SWIG
  public ATransformWithAutoDiff<TuckeyGH>
#else
  public ATransform
#endif
{
public:
  TuckeyGH(double g = 0.0, double h = 0.0);
  TuckeyGH(const TuckeyGH& r)            = default;
  TuckeyGH& operator=(const TuckeyGH& r) = default;
  virtual ~TuckeyGH()                    = default;
  IMPLEMENT_CLONING(TuckeyGH)
  bool hasParameters() const override { return true; }
  double getGValue() const { return _G.getValue(); }
  void setGValue(double g) { _G.setValue(g); }
  double getHValue() const { return _H.getValue(); }
  void setHValue(double h) { _H.setValue(h); }
#ifndef SWIG
  void initParams(double min = 0, double max = INF) override;
#endif
  void appendParams(ListParams& listParams) override;

  String getName() const override { return "TuckeyGH"; }

  template<typename T>
  T evalImpl(T h) const
  {
    double g  = _G.getValue();
    double hv = _H.getValue();
    if (g == 0.)
    {
      return h * exp(0.5 * hv * h * h);
    }

    return (exp(g * h) - 1.) / g * exp(0.5 * hv * h * h);
  }

  VectorDouble getParams() const override;
private:
  ParamInfo _G;
  ParamInfo _H;
};

} // namespace gstlrn
