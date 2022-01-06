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
#include "Basic/AStringFormat.hpp"

AStringFormat::AStringFormat(int level)
    : _level(level)
{
}

AStringFormat::AStringFormat(const AStringFormat& r)
    : _level(r._level)
{
}

AStringFormat& AStringFormat::operator=(const AStringFormat& r)
{
  if (this != &r)
  {
    _level = r._level;
  }
  return *this;
}


AStringFormat::~AStringFormat()
{
}

