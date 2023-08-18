/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Basic/FunctionalSpirale.hpp"
#include "Basic/AFunctional.hpp"
#define _USE_MATH_DEFINES // To make M_PI available under windows
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

FunctionalSpirale::FunctionalSpirale(double a,
                                     double b,
                                     double c,
                                     double d,
                                     double sx,
                                     double sy)
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

/**
 * return the angle of the spiral at a given coordinate
 * @param coor 2-D coordinates of the target
 * @return
 */
double FunctionalSpirale::getFunctionValue(const VectorDouble& coor) const
{
  double x = coor[0] - _sx;
  double y = coor[1] - _sy;
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

/**
 * return the anisotropy rotation matrix at a given coordinate
 * @param coor 2-D coordinates of the target
 * @return
 */
VectorVectorDouble FunctionalSpirale::getFunctionVectors(const VectorDouble& coor) const
{
  double x = coor[0] - _sx;
  double y = coor[1] - _sy;
  double u1 = _linearCombination(x, y, _a, _b);
  double u2 = _linearCombination(x, y, _c, _d);
  double norm = sqrt(u1 * u1 + u2 * u2);
  u1 /= norm;
  u2 /= norm;
  VectorDouble vec1 = { u1,u2};
  VectorDouble vec2 = {-u2,u1};

  VectorVectorDouble vec = { vec1, vec2};
  return vec;
}

