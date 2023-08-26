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
#include "LithoRule/RuleStringFormat.hpp"

RuleStringFormat::RuleStringFormat(int level)
    : AStringFormat(level),
      _flagProp(false),
      _flagThresh(false)
{
  if (level > 1)
  {
    _flagProp = true;
    _flagThresh = true;
  }
}

RuleStringFormat::RuleStringFormat(const RuleStringFormat& r)
    : AStringFormat(r),
      _flagProp(r._flagProp),
      _flagThresh(r._flagThresh)
{
}

RuleStringFormat& RuleStringFormat::operator=(const RuleStringFormat& r)
{
  if (this != &r)
  {
    AStringFormat::operator=(r);
    _flagProp = r._flagProp;
    _flagThresh = r._flagThresh;
  }
  return *this;
}

RuleStringFormat::~RuleStringFormat()
{
}
