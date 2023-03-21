/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
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

