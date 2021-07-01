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
#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"

ADriftElem::ADriftElem(const ENUM_DRIFTS& type,
                       const CovContext& ctxt,
                       int rankFex)
    : ADrift(ctxt.getSpace()),
      _ctxt(ctxt),
      _type(type),
      _rankFex(rankFex),
      _orderIRF(0)
{
}

ADriftElem::ADriftElem(const ADriftElem &r)
    : ADrift(r),
      _ctxt(r._ctxt),
      _type(r._type),
      _rankFex(r._rankFex),
      _orderIRF(r._orderIRF)
{
}

ADriftElem& ADriftElem::operator=(const ADriftElem &r)
{
  if (this != &r)
  {
    ADrift::operator=(r);
    _ctxt = r._ctxt;
    _type = r._type;
    _rankFex = r._rankFex;
    _orderIRF = r._orderIRF;
  }
  return *this;
}

ADriftElem::~ADriftElem()
{
}

bool ADriftElem::isConsistent(const ASpace* space) const
{
  return true;
}

std::string ADriftElem::toString(int level) const
{
  std::stringstream sstr;
  sstr << getDriftName();
  if (getType() == DRIFT_F)
    sstr << " - Rank=" << getRankFex();
  return sstr.str();
}
