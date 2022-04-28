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
#include "Boolean/AToken.hpp"
#include "Basic/AException.hpp"

AToken::AToken()
    : AStringable(),
      _factorX2Y(1.),
      _factorX2Z(1.),
      _factorY2Z(1.),
      _proportion(1.),
      _paramNames(),
      _params()
{
}

AToken::AToken(const AToken &r)
    : AStringable(r),
      _factorX2Y(r._factorX2Y),
      _factorX2Z(r._factorX2Z),
      _factorY2Z(r._factorY2Z),
      _proportion(r._proportion),
      _paramNames(r._paramNames),
      _params(r._params)
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
    _paramNames = r._paramNames;
    _params = r._params;
  }
  return *this;
}

AToken::~AToken()
{
}

String AToken::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "Token:"<< getType().getDescr() << " - Nb. params=" << getNParams() <<
      " - Proportion=" <<_proportion << std::endl;

  for (int ipar = 0; ipar < getNParams(); ipar++)
  {
    sstr << getParamName(ipar) << getParam(ipar).toString() << std::endl;
  }

  return sstr.str();
}

void AToken::initParams(int count)
{
  _paramNames.resize(count);
  _params.resize(count);

  for (int ipar = 0; ipar < count; ipar++)
  {
    _params[ipar] = TokenParameter();
  }
}

void AToken::setLaw(int ipar, ETLaw law)
{
  if (! _isValidParamIndex(ipar)) return;
  _params[ipar].setLaw(law);
}

void AToken::setValarg(int ipar, int iarg, double value)
{
  if (! _isValidParamIndex(ipar)) return;
  _params[ipar].setValarg(iarg, value);
}

void AToken::setParamName(int ipar, const String& name)
{
  if (! _isValidParamIndex(ipar)) return;
  _paramNames[ipar] = name;
}

String AToken::getParamName(int ipar) const
{
  if (! _isValidParamIndex(ipar)) return String();
  return _paramNames[ipar];
}

const TokenParameter& AToken::getParam(int ipar) const
{
  if (! _isValidParamIndex(ipar))
    my_throw("Argument invalid");
  return _params[ipar];
}

double AToken::generateParam(int ipar) const
{
  if (! _isValidParamIndex(ipar)) return TEST;
 return _params[ipar].generateValue();
}

bool AToken::_isValidParamIndex(int ipar) const
{
  int npar = (int) _params.size();
  if (ipar < 0 || ipar >= npar)
  {
    messerr("fdsfsd");
    return false;
  }
  return true;
}

