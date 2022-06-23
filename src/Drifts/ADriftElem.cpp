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
#include "Drifts/ADriftElem.hpp"
#include "Drifts/DriftFactory.hpp"
#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Drifts/EDrift.hpp"

ADriftElem::ADriftElem(const EDrift& type,
                       const CovContext& ctxt,
                       int rankFex)
    : ADrift(ctxt.getSpace()), /// TODO : shared pointer
      ASerializable(),
      _ctxt(ctxt),
      _type(type),
      _rankFex(rankFex)
{
}

ADriftElem::ADriftElem(const ADriftElem &r)
    : ADrift(r),
      ASerializable(r),
      _ctxt(r._ctxt), /// TODO : shared pointer
      _type(r._type),
      _rankFex(r._rankFex)
{
}

ADriftElem& ADriftElem::operator=(const ADriftElem &r)
{
  if (this != &r)
  {
    ADrift::operator=(r);
    ASerializable::operator=(r);
    _ctxt = r._ctxt;
    _type = r._type;
    _rankFex = r._rankFex;
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
  if (getType() == EDrift::F)
    sstr << " - Rank=" << getRankFex();
  return sstr.str();
}

bool ADriftElem::_deserialize(std::istream& is, bool /*verbose*/)
{
  bool ret = true;
  int type;
  ret = ret && _recordRead<int>(is, "Drift Function", type);
  _type = EDrift::fromValue(type);
  _rankFex = 0;
  return ret;
}

bool ADriftElem::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os,"Drift characteristics", getType().getValue());
  return ret;
}

