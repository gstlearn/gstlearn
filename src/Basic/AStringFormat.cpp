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
#include "Basic/AStringFormat.hpp"

AStringFormat::AStringFormat(int level)
  : _level(level)
  , _title()
{
}

AStringFormat::AStringFormat(const String& title)
  : _level()
  , _title(title)
{
}

AStringFormat::AStringFormat(const AStringFormat& r)
  : _level(r._level)
  , _title(r._title)
{
}

AStringFormat& AStringFormat::operator=(const AStringFormat& r)
{
  if (this != &r)
  {
    _level = r._level;
    _title = r._title;
  }
  return *this;
}

AStringFormat::~AStringFormat()
{
}

