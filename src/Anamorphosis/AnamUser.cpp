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
#include "Anamorphosis/AnamUser.hpp"

#include "geoslib_f.h"

#include "Basic/AException.hpp"

AnamUser::AnamUser()
    : AnamContinuous(),
      _y2z_function(nullptr),
      _z2y_function(nullptr)
{
}

AnamUser::AnamUser(const AnamUser &m)
    : AnamContinuous(m),
      _y2z_function(m._y2z_function),
      _z2y_function(m._z2y_function)
{

}

AnamUser& AnamUser::operator=(const AnamUser &m)
{
  if (this != &m)
  {
    _y2z_function = m._y2z_function;
    _z2y_function = m._z2y_function;
  }
  return *this;
}

AnamUser::~AnamUser()
{

}

String AnamUser::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "User defined Anamorphosis" << std::endl;

  return sstr.str();
}

void AnamUser::calculateMeanAndVariance()
{
  messerr("AnamUser: This funtion does not make sense");
}

int AnamUser::_deserialize(std::istream& /*is*/, bool /*verbose*/)
{
  messerr("AnamUser: Cannot be deserialized");
  return 1;
}

int AnamUser::_serialize(std::ostream& /*os*/, bool /*verbose*/) const
{
  messerr("AnamUser: Cannot be serialized");
  return 1;
}

double AnamUser::GaussianToRawValue(double h) const
{
  if (_y2z_function == nullptr) return TEST;
  return _y2z_function(h);
}

double AnamUser::RawToGaussianValue(double h) const
{
  if (_z2y_function == nullptr) return TEST;
  return _z2y_function(h);
}
