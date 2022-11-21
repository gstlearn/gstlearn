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
#include "Basic/AStringable.hpp"
#include "Simulation/SimuFFTParam.hpp"

#include <math.h>

SimuFFTParam::SimuFFTParam(bool flag_aliasing,
                           double percent)
    : AStringable(),
      _flagAliasing(flag_aliasing),
      _percent(percent)
{
}

SimuFFTParam::SimuFFTParam(const SimuFFTParam &r)
    : AStringable(r),
      _flagAliasing(r._flagAliasing),
      _percent(r._percent)
{
}

SimuFFTParam& SimuFFTParam::operator=(const SimuFFTParam &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
    _flagAliasing = r._flagAliasing;
    _percent = r._percent;
  }
  return *this;
}

SimuFFTParam::~SimuFFTParam()
{
}

String SimuFFTParam::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  if (_flagAliasing)
    sstr << "Perform intermediate mesh discretization in order to reduce aliasing" << std::endl;
  sstr << "Percentage of Covariance used for field extension" << std::endl;

  return sstr.str();
}
