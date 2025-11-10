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

namespace gstlrn
{

TuckeyGH::TuckeyGH(double g, double h)
  : ATransformWithAutoDiff<TuckeyGH>()
  , _G(ParamInfo("Tuckey G-Value", g, {-INF, INF}))
  , _H(ParamInfo("Tuckey H-Value", h, {.0, 1})) // H > 0 guarantees monotonicity
{
}

void TuckeyGH::initParams()
{
  _G.setValue(0.);
  //_H.setValue(EPSILON4); // H not on its lower bound 
}


void TuckeyGH::appendParams(ListParams& listParams)
{
  listParams.addParam(_G);
  //listParams.addParam(_H);
}

} // namespace gstlrn
