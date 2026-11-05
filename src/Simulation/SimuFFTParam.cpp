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
#include "Simulation/SimuFFTParam.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Message.hpp"

#include <cmath>

namespace gstlrn
{
  SimuFFTParam::SimuFFTParam(bool flag_aliasing, double percent)
    : AStringable()
    , _flagAliasing()
    , _percent(percent)
  {
    _checkAliasing(flag_aliasing);
  }

  SimuFFTParam::SimuFFTParam(const SimuFFTParam& r)
    : AStringable(r)
    , _flagAliasing(r._flagAliasing)
    , _percent(r._percent)
  {
  }

  SimuFFTParam& SimuFFTParam::operator=(const SimuFFTParam& r)
  {
    if (this != &r)
    {
      AStringable::operator=(r);
      _flagAliasing = r._flagAliasing;
      _percent = r._percent;
    }
    return *this;
  }

  SimuFFTParam::~SimuFFTParam() {}

  String SimuFFTParam::toString(const AStringFormat* /*strfmt*/) const
  {
    std::stringstream sstr;

    if (_flagAliasing)
      sstr << "Perform intermediate mesh discretization in order to reduce "
              "aliasing (OBSOLETE)"
           << std::endl;
    sstr << "Percentage of Covariance used for field extension (1% minimum)"
         << std::endl;

    return sstr.str();
  }

  void SimuFFTParam::_checkAliasing(bool aliasing)
  {
    if (!aliasing) return;
    messerr(
      "The 'aliasing' feature is now OBSOLETE for the Simulation by FFT "
      "method");
  }
} // namespace gstlrn
