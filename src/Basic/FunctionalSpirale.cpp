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
#include "Basic/FunctionalSpirale.hpp"
#include "Basic/AFunctional.hpp"
#include <math.h>

FunctionalSpirale::FunctionalSpirale()
    : AFunctional(2),
      _a(0.),
      _b(0.),
      _c(0.),
      _d(0.),
      _sx(0.),
      _sy(0.)
{
}

FunctionalSpirale::FunctionalSpirale(double a, double b, double c, double d, double sx, double sy)
    : AFunctional(2),
      _a(a),
      _b(b),
      _c(c),
      _d(d),
      _sx(sx),
      _sy(sy)
{
}

FunctionalSpirale::FunctionalSpirale(const FunctionalSpirale &m)
    : AFunctional(m),
      _a(m._a),
      _b(m._b),
      _c(m._c),
      _d(m._d),
      _sx(m._sx),
      _sy(m._sy)
{
}

FunctionalSpirale& FunctionalSpirale::operator=(const FunctionalSpirale &m)
{
  if (this != &m)
  {
    AFunctional::operator=(m);
    _a = m._a;
    _b = m._b;
    _c = m._c;
    _d = m._d;
    _sx = m._sx;
    _sy = m._sy;
  }
  return *this;
}

FunctionalSpirale::~FunctionalSpirale()
{
}

double FunctionalSpirale::_linearCombination(double x, double y, double a, double b) const
{
    return a*x + b*y;
}

double FunctionalSpirale::getFunctionValue(const VectorDouble& pos) const
{
  double x = pos[0] - _sx;
  double y = pos[1] - _sy;
  double u1 = _linearCombination(x, y, _a, _b);
  double u2 = _linearCombination(x, y, _c, _d);
  double norm = sqrt(u1 * u1 + u2 * u2);
  if (norm > 0)
  {
    return acos(u2 / norm) / M_PI * 180. * ((u1 >= 0) ? 1. : -1.);
  }
  else
  {
    return 0.;
  }
}

