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
#include "../../include/Boolean/AToken.hpp"

AToken::AToken()
    : AStringable(),
      _factorX2Y(1.),
      _factorX2Z(1.),
      _factorY2Z(1.),
      _proportion(1.)
{
}

AToken::AToken(const AToken &r)
    : AStringable(r),
      _factorX2Y(r._factorX2Y),
      _factorX2Z(r._factorX2Z),
      _factorY2Z(r._factorY2Z),
      _proportion(r._proportion)
{
}

AToken& AToken::operator=(const AToken &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
    _factorX2Y = r._factorX2Y;
    _factorX2Z = r._factorX2Z;
    _factorY2Z = r._factorY2Z;
    _proportion = r._proportion;
  }
  return *this;
}

AToken::~AToken()
{
}

String AToken::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "Token:"<< getType().getDescr() << " - Nb. params=" << getNArgs() <<
      " - Proportion=" <<_proportion << std::endl;

  for (int i = 0; i < getNArgs(); i++)
  {
    sstr << getParamName(i) << getParam(i).toString() << std::endl;
  }

  return sstr.str();
}
