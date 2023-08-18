/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include <Geometry/BiTargetCheckCode.hpp>
#include "geoslib_f.h"

#include "Space/SpaceTarget.hpp"

BiTargetCheckCode::BiTargetCheckCode(int optcode, double tolcode)
    : ABiTargetCheck(),
      _optCode(optcode),
      _tolCode(tolcode)
{
}

BiTargetCheckCode::BiTargetCheckCode(const BiTargetCheckCode &r)
    : ABiTargetCheck(r),
      _optCode(r._optCode),
      _tolCode(r._tolCode)
{
}

BiTargetCheckCode& BiTargetCheckCode::operator=(const BiTargetCheckCode &r)
{
  if (this != &r)
  {
    ABiTargetCheck::operator=(r);
    _optCode = r._optCode;
    _tolCode = r._tolCode;
  }
  return *this;
}

BiTargetCheckCode::~BiTargetCheckCode()
{
}

String BiTargetCheckCode::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  if (_optCode == 1)
    sstr << "- Must have similar Codes (tol=" << _tolCode << ")" << std::endl;
  else
    sstr << "- Must have different Codes" << std::endl;

  return sstr.str();
}

bool BiTargetCheckCode::isOK(const SpaceTarget &T1, const SpaceTarget &T2) const
{
  double code1 = T1.getCode();
  double code2 = T2.getCode();

  switch (_optCode)
  {
    case 1: /* Code must be close */
      if (ABS(code1 - code2) > _tolCode) return false;
      break;

    case 2: /* Code must be different */
      if (code1 == code2) return false;
      break;
  }
  return true;
}
