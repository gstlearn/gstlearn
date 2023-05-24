/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"

#include "Geometry/ABiPointCheck.hpp"

ABiPointCheck::ABiPointCheck()
    : AStringable()
{
}

ABiPointCheck::ABiPointCheck(const ABiPointCheck &r)
    : AStringable(r)
{
}

ABiPointCheck& ABiPointCheck::operator=(const ABiPointCheck &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
  }
  return *this;
}

ABiPointCheck::~ABiPointCheck()
{
}

String ABiPointCheck::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "BiPoint Check" << std::endl;

  return sstr.str();
}

