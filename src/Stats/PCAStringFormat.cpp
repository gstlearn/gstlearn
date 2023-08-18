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
#include "Stats/PCAStringFormat.hpp"

PCAStringFormat::PCAStringFormat(int level)
    : AStringFormat(level),
      _flagCenter(true),
      _flagStats(true)
{
  if (level == 0)
  {
    _flagCenter = false;
    _flagStats = false;
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
