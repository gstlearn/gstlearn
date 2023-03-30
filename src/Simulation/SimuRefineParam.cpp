/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
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

