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

#include "Drifts/ADriftElem.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"

ADriftElem::ADriftElem(const CovContext& ctxt)
    : ADrift(ctxt.getSpace()), /// TODO : shared pointer
      _ctxt(ctxt)
{
}

ADriftElem::ADriftElem(const ADriftElem &r)
    : ADrift(r),
      _ctxt(r._ctxt) /// TODO : shared pointer
{
}

ADriftElem& ADriftElem::operator=(const ADriftElem &r)
{
  if (this != &r)
  {
    ADrift::operator=(r);
    _ctxt = r._ctxt;
  }
  return *this;
}

ADriftElem::~ADriftElem()
{
}

bool ADriftElem::isConsistent(const ASpace* /*space*/) const
{
  return true;
}

String ADriftElem::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << getDriftName();
  return sstr.str();
}
