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
#include "Anamorphosis/AnamUser.hpp"

namespace gstlrn
{
AnamUser::AnamUser()
  : AnamContinuous()
  , _y2z_function(nullptr)
  , _z2y_function(nullptr)
{
}

AnamUser::AnamUser(const AnamUser& m)
  : AnamContinuous(m)
  , _y2z_function(m._y2z_function)
  , _z2y_function(m._z2y_function)
{
}

AnamUser& AnamUser::operator=(const AnamUser& m)
{
  if (this != &m)
  {
    AnamContinuous::operator=(m);
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

bool AnamUser::_deserializeAscii(std::istream& /*is*/)
{
  messerr("AnamUser: Cannot be deserialized");
  return false;
}

bool AnamUser::_serializeAscii(std::ostream& /*os*/) const
{
  messerr("AnamUser: Cannot be serialized");
  return false;
}

double AnamUser::transformToRawValue(double h) const
{
  if (_y2z_function == nullptr) return TEST;
  return _y2z_function(h);
}

double AnamUser::rawToTransformValue(double h) const
{
  if (_z2y_function == nullptr) return TEST;
  return _z2y_function(h);
}
#ifdef HDF5
bool AnamUser::deserializeH5(H5::Group& grp)
{
  DECLARE_UNUSED(grp);
  messerr("AnamUser: Cannot be deserialized");
  return false;
}

bool AnamUser::serializeH5(H5::Group& grp) const
{
  DECLARE_UNUSED(grp);
  messerr("AnamUser: Cannot be serialized");
  return false;
}
#endif
} // namespace gstlrn