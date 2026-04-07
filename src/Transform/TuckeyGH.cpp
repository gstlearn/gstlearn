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

#include "Transform/TuckeyGH.hpp"
#include "Transform/ATransform.hpp"
#include "geoslib_define.h"

namespace gstlrn
{

  TuckeyGH::TuckeyGH(double g, double h)
    : ATransformWithAutoDiff<TuckeyGH>()
    , _G(ParamInfo("Tuckey G-Value", g, {-INF, INF}))
    , _H(ParamInfo(
        "Tuckey H-Value",
        h,
        {.0, 1})) // H > 0 guarantees monotonicity
  {
  }

  void TuckeyGH::initParams(double min, double max)
  {
    DECLARE_UNUSED(min, max)
    _G.setValue(0.);
    _H.setValue(EPSILON4); // H not on its lower bound
  }

  void TuckeyGH::appendParams(ListParams& listParams)
  {
    listParams.addParam(_G);
    // listParams.addParam(_H);
  }

  VectorDouble TuckeyGH::getParams() const
  {
    return {_G.getValue(), _H.getValue()};
  } // namespace gstlrn

} // namespace gstlrn
