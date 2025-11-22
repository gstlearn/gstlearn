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
#include "Basic/AStringable.hpp"

#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <typeinfo>

namespace gstlrn
{
AStringable::AStringable()
{
}

/**
 * Copy constructor: don't copy temporary file info
 */
AStringable::AStringable(const AStringable& /*r*/)
{
}
/**
 * Assignment operator: don't copy temporary file info
 */
AStringable& AStringable::operator=(const AStringable& /*r*/)
{
  return *this;
}

AStringable::~AStringable()
{
}

String AStringable::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << "toString is not yet implemented for " << typeid(*this).name() << std::endl;
  return sstr.str();
}

/**
 * Send the String to the display function
 */
void AStringable::display(const AStringFormat* strfmt) const
{
  if (strfmt != nullptr)
  {
    if (strfmt->hasTitle())
    {
      message_extern(strfmt->getTitle().c_str());
      message_extern("\n");
    }
  }
  message_extern(toString(strfmt).c_str());
}

void AStringable::display(Id level) const
{
  AStringFormat sf(level);
  display(&sf);
}

} // namespace gstlrn
