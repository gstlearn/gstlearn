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
#include "geoslib_f_private.h"

#include "Drifts/ADrift.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"

ADrift::ADrift()
    : AStringable()
{
}

ADrift::ADrift(const ADrift &r)
    : AStringable(r)
{
}

ADrift& ADrift::operator=(const ADrift &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
  }
  return *this;
}

ADrift::~ADrift()
{
}

String ADrift::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << getDriftName();
  return sstr.str();
}
