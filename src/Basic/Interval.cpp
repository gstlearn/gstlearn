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
#include "Basic/Interval.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

Interval::Interval(double vmin, double vmax, bool mininc, bool maxinc)
    : _vmin(vmin),
      _vmax(vmax),
      _minIncluded(mininc),
      _maxIncluded(maxinc)
{
  _modifyUnbounded();
  if (! _isValid()) throw("Interval is not valid");
}

Interval::Interval(const Interval &m)
    : _vmin(m._vmin),
      _vmax(m._vmax),
      _minIncluded(m._minIncluded),
      _maxIncluded(m._maxIncluded)
{

}

Interval& Interval::operator=(const Interval &m)
{
  if (this != &m)
  {
    _vmin = m._vmin;
    _vmax = m._vmax;
    _minIncluded = m._minIncluded;
    _maxIncluded = m._maxIncluded;
  }
  return *this;
}

Interval::~Interval()
{

}

String Interval::toString(int level) const
{
  std::stringstream sstr;
  if (FFFF(_vmin))
    sstr << "] -Inf";
  else
  {
    if (_minIncluded)
      sstr << "[ " << _vmin;
    else
      sstr << " ]" << _vmin;
  }
  sstr << " ; ";
  if (FFFF(_vmax))
    sstr << " +Inf [";
  else
  {
    if (_maxIncluded)
      sstr << _vmax << " ]";
    else
      sstr << _vmax << " [";
  }
  return sstr.str();
}

void Interval::init(double vmin, double vmax, bool mininc, bool maxinc)
{
  _vmin = vmin;
  _vmax = vmax;
  _minIncluded = mininc;
  _maxIncluded = maxinc;
  _modifyUnbounded();
  if (! _isValid()) throw("Interval is not valid");
}

bool Interval::isInside(double value) const
{
  if (FFFF(value)) return false;
  if (isBelow(value)) return false;
  if (isAbove(value)) return false;
  return true;
}

bool Interval::isBelow(double value) const
{
  if (FFFF(value)) return false;
  if (_minIncluded)
  {
    if (! FFFF(_vmin) && value > _vmin) return false;
  }
  else
  {
    if (! FFFF(_vmin) && value >= _vmin) return false;
  }
  return true;
}

bool Interval::isAbove(double value) const
{
  if (FFFF(value)) return false;
  if (_maxIncluded)
  {
    if (! FFFF(_vmax) && value < _vmax) return false;
  }
  else
  {
    if (! FFFF(_vmax) && value <= _vmax) return false;
  }
  return true;
}

VectorDouble Interval::getBounds() const
{
  VectorDouble vec(2);
  vec[0] = _vmin;
  vec[1] = _vmax;
  return vec;
}

bool Interval::isValid() const
{
  if (_vmin < _vmax) return true;
  if (_vmin > _vmax) return false;
  if (_minIncluded && _maxIncluded) return true;
  return false;
}

bool Interval::isDisjoint(const Interval& m) const
{
  if (m._vmin < _vmax)
  {
    // The Interval 'm' is located below 'this'

    if (m._vmax < _vmin) return true;
    if (m._vmax > _vmin) return false;

    if (m._maxIncluded && ! _minIncluded) return true;
    if (! m._maxIncluded && _minIncluded) return true;
    return false;
  }
  else
  {
    // The Interval 'm' is located above 'this'

    if (m._minIncluded && ! _maxIncluded) return true;
    if (! m._minIncluded && _maxIncluded) return true;
    return false;
  }
  return false;
}

void Interval::_modifyUnbounded(void)
{
  if (FFFF(_vmin)) _minIncluded = false;
  if (FFFF(_vmax)) _maxIncluded = false;
}

bool Interval::_isValid(void)
{
  if (FFFF(_vmin) || FFFF(_vmax)) return true;
  if (_vmin < _vmax) return true;
  if (_vmin > _vmax)
  {
    messerr("Interval Definition: Lower Bound(%lf) should be smaller than Upper Bound(%lf)",
            _vmin,_vmax);
    return false;
  }

  if (! _minIncluded ||  ! _maxIncluded)
  {
    messerr("Interval Definition: Bounds are equal; then interval should be closed");
    return false;
  }
  return true;
}
