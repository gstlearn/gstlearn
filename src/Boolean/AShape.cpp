/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Boolean/AShape.hpp"
#include "Basic/AException.hpp"

AShape::AShape()
    : AStringable(),
      _factorX2Y(0.),
      _factorX2Z(0.),
      _factorY2Z(0.),
      _proportion(1.),
      _paramNames(),
      _params()
{
}

AShape::AShape(const AShape &r)
    : AStringable(r),
      _factorX2Y(r._factorX2Y),
      _factorX2Z(r._factorX2Z),
      _factorY2Z(r._factorY2Z),
      _proportion(r._proportion),
      _paramNames(r._paramNames),
      _params(r._params)
{
}

AShape& AShape::operator=(const AShape &r)
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

AShape::~AShape()
{
}

String AShape::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << getType().getDescr() << " - Proportion=" << _proportion
      << std::endl;

  for (int ipar = 0; ipar < getNParams(); ipar++)
  {
    sstr << "- " << getParamName(ipar) << ":" << getParam(ipar).toString();
  }

  if (_factorX2Y > 0.)
    sstr << "Y-Extension = X_Extension * "<< _factorX2Y << std::endl;
  if (_factorX2Z > 0.)
    sstr << "Z-Extension = X_Extension * "<< _factorX2Z << std::endl;
  if (_factorY2Z > 0.)
    sstr << "Z-Extension = Y_Extension * "<< _factorY2Z << std::endl;

  return sstr.str();
}

void AShape::initParams(int count)
{
  _paramNames.resize(count);
  _params.resize(count);

  for (int ipar = 0; ipar < count; ipar++)
  {
    _params[ipar] = ShapeParameter();
  }
}

void AShape::setLaw(int ipar, ELaw law)
{
  if (! _isValidParamIndex(ipar)) return;
  _params[ipar].setLaw(law);
}

void AShape::setParam(int ipar, int iarg, double value)
{
  if (! _isValidParamIndex(ipar)) return;
  _params[ipar].setValarg(iarg, value);
}

void AShape::setParamName(int ipar, const String& name)
{
  if (! _isValidParamIndex(ipar)) return;
  _paramNames[ipar] = name;
}

void AShape::setParamDefault(int ipar,
                             const String& name,
                             double value)
{
  if (! _isValidParamIndex(ipar)) return;
  _paramNames[ipar] = name;
  _params[ipar].setValarg(0, value);
}

String AShape::getParamName(int ipar) const
{
  if (! _isValidParamIndex(ipar)) return String();
  return _paramNames[ipar];
}

double AShape::getParam(int ipar, int iarg) const
{
  if (! _isValidParamIndex(ipar)) return TEST;
  return _params[ipar].getValarg(iarg);
}

const ShapeParameter& AShape::getParam(int ipar) const
{
  if (! _isValidParamIndex(ipar))
    my_throw("Argument invalid");
  return _params[ipar];
}

double AShape::generateParam(int ipar) const
{
  if (! _isValidParamIndex(ipar)) return TEST;
 return _params[ipar].generateValue();
}

bool AShape::_isValidParamIndex(int ipar) const
{
  int npar = (int) _params.size();
  if (ipar < 0 || ipar >= npar)
  {
    messerr("Index %d is not valid. It should lie in [0,%d[",
            ipar,npar);
    return false;
  }
  return true;
}

