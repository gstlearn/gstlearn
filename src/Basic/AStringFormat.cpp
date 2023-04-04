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

