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
#include "Stats/PCAStringFormat.hpp"

PCAStringFormat::PCAStringFormat(int level)
    : AStringFormat(level),
      _flagCenter(false),
      _flagStats(false)
{
  if (level > 1)
  {
    _flagCenter = true;
    _flagStats = true;
  }
}

PCAStringFormat::PCAStringFormat(const PCAStringFormat& r)
    : AStringFormat(r),
      _flagCenter(r._flagCenter),
      _flagStats(r._flagStats)
{
}

PCAStringFormat& PCAStringFormat::operator=(const PCAStringFormat& r)
{
  if (this != &r)
  {
    AStringFormat::operator=(r);
    _flagCenter = r._flagCenter;
    _flagStats = r._flagStats;
  }
  return *this;
}

PCAStringFormat::~PCAStringFormat()
{
}
