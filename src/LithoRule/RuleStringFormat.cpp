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
