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
#include "geoslib_enum.h"

#include "Anamorphosis/PPMT.hpp"
#include "Db/Db.hpp"
#include "Basic/AException.hpp"

PPMT::PPMT()
    : AStringable(),
      _anams(),
      _directions()
{
}

PPMT::PPMT(const PPMT &m)
    : AStringable(m),
      _anams(m._anams),
      _directions(m._directions)
{
}

PPMT& PPMT::operator=(const PPMT &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _anams = m._anams;
    _directions = m._directions;
  }
  return *this;
}

PPMT::~PPMT()
{
}

String PPMT::toString(const AStringFormat* strfmt) const
{
  SYMBOL_UNUSED(strfmt);

  std::stringstream sstr;

  int niter = (int) _anams.size();
  for (int iter = 0; iter < niter; iter++)
  {
    sstr << _anams[iter].toString(strfmt);
    sstr << "Direction = " << _directions[iter] << std::endl;
  }

  return sstr.str();
}

void PPMT::addIteration(const AnamHermite& anam, const VectorDouble& dir)
{
  _anams.push_back(anam);
  _directions.push_back(dir);
}
