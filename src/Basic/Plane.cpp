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
#include "Basic/Plane.hpp"

Plane::Plane()
    : AStringable(),
      _coor(3),
      _intercept(0),
      _value(0.),
      _rndval(0.)
{
}

Plane::Plane(const Plane &m)
    : AStringable(m),
      _coor(m._coor),
      _intercept(m._intercept),
      _value(m._value),
      _rndval(m._rndval)
{
}

Plane& Plane::operator=(const Plane &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _coor = m._coor;
    _intercept = m._intercept;
    _value = m._value;
    _rndval = m._rndval;
  }
  return *this;
}

Plane::~Plane()
{

}

String Plane::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;


  return sstr.str();
}
