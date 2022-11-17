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
#include "Simulation/SimuRefineParam.hpp"

#include <math.h>

SimuRefineParam::SimuRefineParam(int nmult, bool flag_SK)
    : AStringable(),
      _nmult(nmult),
      _flagSK(flag_SK)
{
}

SimuRefineParam::SimuRefineParam(const SimuRefineParam &r)
    : AStringable(r),
      _nmult(r._nmult),
      _flagSK(r._flagSK)
{
}

SimuRefineParam& SimuRefineParam::operator=(const SimuRefineParam &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
    _nmult = r._nmult;
    _flagSK = r._flagSK ;
  }
  return *this;
}

SimuRefineParam::~SimuRefineParam()
{
}

String SimuRefineParam::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "Refinement Factor = " << _nmult << std::endl;
  if (_flagSK)
    sstr << "Refinement using Simple Kriging" << std::endl;
  else
    sstr << "Refinement using Ordinary Kriging" << std::endl;

  return sstr.str();
}

