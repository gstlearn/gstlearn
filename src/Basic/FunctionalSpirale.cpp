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
#include "geoslib_define.h"

#include "Db/Db.hpp"
#include "Covariances/CovAniso.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "Basic/AFunctional.hpp"
#include <math.h>

FunctionalSpirale::FunctionalSpirale()
    : AFunctional(2),
      _a(0.),
      _b(0.),
      _c(0.),
      _d(0.),
      _xcenter(0.),
      _ycenter(0.)
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
      _xcenter(sx),
      _ycenter(sy)
{
}

FunctionalSpirale::FunctionalSpirale(const FunctionalSpirale &m)
    : AFunctional(m),
      _a(m._a),
      _b(m._b),
      _c(m._c),
      _d(m._d),
      _xcenter(m._xcenter),
      _ycenter(m._ycenter)
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
    _xcenter = m._xcenter;
    _ycenter = m._ycenter;
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
 * Return the angle (in degrees) of the spiral at a given coordinate
 * @param coor 2-D coordinates of the target
 * @return
 */
double FunctionalSpirale::getFunctionValue(const VectorDouble& coor) const
{
  double x = coor[0] - _xcenter;
  double y = coor[1] - _ycenter;
  double u1 = _linearCombination(x, y, _a, _b);
  double u2 = _linearCombination(x, y, _c, _d);
  double norm = sqrt(u1 * u1 + u2 * u2);
  if (norm > 0)
  {
    double a2ndeg = acos(u2 / norm) * 180. / GV_PI;
    return (u1 >= 0) ? -a2ndeg : a2ndeg;
  }
  else
  {
    return 0.;
  }
}

/**
 * Return the anisotropy rotation matrix at a given coordinate
 * @param coor 2-D coordinates of the target
 * @return
 */
MatrixSquareGeneral FunctionalSpirale::getFunctionMatrix(const VectorDouble& coor) const
{
  int ndim = 2;
  MatrixSquareGeneral dirs = MatrixSquareGeneral(ndim);

  double angle = getFunctionValue(coor) * GV_PI / 180.;
  double u1 = cos(angle);
  double u2 = sin(angle);
  dirs.setValue(0, 0,  u1);
  dirs.setValue(1, 0, -u2);
  dirs.setValue(0, 1,  u2);
  dirs.setValue(1, 1,  u1);
  return dirs;
}

VectorVectorDouble FunctionalSpirale::getFunctionVectors(const Db *db, const CovAniso* cova) const
{
  if (db == nullptr) return VectorVectorDouble();
  if (getNdim() != db->getNDim())
  {
    messerr("You cannot evaluate the function on input Db: they do not have the same Space Dimension");
    return VectorVectorDouble();
  }

  int nech = db->getSampleNumber();
  VectorVectorDouble vec(3);
  vec[0].resize(nech);
  vec[1].resize(nech);
  vec[2].resize(nech);

  MatrixSquareSymmetric temp(2);
  MatrixSquareSymmetric hh(2);
  VectorDouble diag = VH::power(cova->getScales(), 2.);
  temp.setDiagonal(diag);

  for (int iech = 0; iech < nech; iech++)
  {
    VectorDouble coor = db->getSampleCoordinates(iech);
    MatrixSquareGeneral rotmat = getFunctionMatrix(coor);
    hh.normMatrix(rotmat, temp);

    vec[0][iech] = hh.getValue(0,0);
    vec[1][iech] = hh.getValue(0,1);
    vec[2][iech] = hh.getValue(1,1);
  }

  return vec;
}
