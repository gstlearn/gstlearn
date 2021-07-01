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
    : AnamContinuous(ANAM_EXTERNAL),
      _y2z_function(nullptr),
      _z2y_function(nullptr)

{
  setType(ANAM_EXTERNAL);
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

String AnamUser::toString(int level) const
{
  std::stringstream sstr;
  sstr << Anam::toString(level);
  sstr << "User defined anamorphosis" << std::endl;
  return sstr.str();
}

void AnamUser::calculateMeanAndVariance()
{
  my_throw("This function is not available for User-defined Anamorphosis");
}
