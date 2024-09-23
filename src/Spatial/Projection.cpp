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
#include "Spatial/Projection.hpp"

#include "Polygon/Polygons.hpp"
#include "Basic/AStringable.hpp"

Projection::Projection(bool flag_mean, double xcenter, double ycenter)
    : _flagMean(flag_mean),
      _xcenter(xcenter),
      _ycenter(ycenter)
{
}

Projection::Projection(bool flag_mean, Db* db)
: _flagMean(flag_mean),
  _xcenter(TEST),
  _ycenter(TEST)
{
  int ndim = db->getLocNumber(ELoc::X);
  if (ndim != 2)
  {
    messerr("The Projection structure is only valid for 2-D space");
    return;
  }

  VectorDouble x = db->getCoordinates(0, true);
  VectorDouble y = db->getCoordinates(1, true);

  _xcenter = VH::mean(x);
  _ycenter = VH::mean(y);
}

Projection::Projection(const Projection &r)
    : _flagMean(r._flagMean),
      _xcenter(r._xcenter),
      _ycenter(r._ycenter)
{
}

Projection& Projection::operator=(const Projection &r)
{
  if (this != &r)
  {
    _flagMean = r._flagMean;
    _xcenter = r._xcenter;
    _ycenter = r._ycenter;
  }
  return *this;
}

Projection::~Projection()
{
}

void Projection::operateInPlace(VectorDouble& coor) const
{
  double xx = coor[0];
  double yy = coor[1];
  if (FFFF(xx) || FFFF(yy)) return;

  double deltax = xx - _xcenter;
  double deltay = yy - _ycenter;

  if (_flagMean)
  {
    // Mean projection

    xx = deltax * 60 * cos(_ycenter * GV_PI / 180.);
    yy = deltay * 60.;
  }
  else
  {
    xx = deltax * 60 * cos((deltay + _ycenter * GV_PI) / 180);
    yy = deltay * 60.;
  }

  coor[0] = xx;
  coor[1] = yy;
}

VectorDouble Projection::operateInvert(const VectorDouble& coor) const
{
  VectorDouble res(2);
  double xx = coor[0];
  double yy = coor[1];
  if (FFFF(xx) || FFFF(yy)) return res;

  if (_flagMean)
  {
    // Mean projection

    xx = _xcenter + xx / (60. * cos(_ycenter * GV_PI / 180.));
    yy = _ycenter + yy / 60.;
  }
  else
  {
    yy = _ycenter + yy / 60.;
    xx = _xcenter + xx / (60. * cos((yy - _ycenter + _ycenter * GV_PI) / 180.));
  }

  res[0] = xx;
  res[1] = yy;
  return res;
}

int Projection::operateVecInPlace(VectorDouble& x, VectorDouble& y) const
{
  int nech = (int) x.size();
  if (nech != (int) y.size())
  {
    messerr("Arguments 'x' and 'y' should have same dimension");
    return 1;
  }

  // Loop on the samples to be projected

  VectorDouble coor(2);
  for (int iech = 0; iech < nech; iech++)
  {
    coor[0] = x[iech];
    coor[1] = y[iech];
    operateInPlace(coor);
    x[iech] = coor[0];
    y[iech] = coor[1];
  }

  return 0;
}

int Projection::operateOnDb(Db *db) const
{
  if (db == nullptr) return 0;
  if (db->getLocNumber(ELoc::X) < 2)
  {
    messerr("This method is dedicated to 2-D space (or more)");
    return 1;
  }

  // Extract the vector of coordinates
  VectorDouble x = db->getCoordinates(0, true);
  VectorDouble y = db->getCoordinates(1, true);

  // Perform the projection
  if (operateVecInPlace(x, y)) return 1;

  // Store the coordinates (in place)
  db->setCoordinates(0, x, true);
  db->setCoordinates(1, y, true);

  return 0;
}

int Projection::operateOnPolygons(Polygons* poly) const
{
  if (poly == nullptr) return 0;
  int npol = poly->getPolyElemNumber();

  // Loop on the Polygon elements
  for (int ipol = 0; ipol < npol; ipol++)
  {
    VectorDouble xx = poly->getX(ipol);
    VectorDouble yy = poly->getY(ipol);

    // Perform the projection
    if (operateVecInPlace(xx, yy)) return 1;

    poly->setX(ipol, xx);
    poly->setY(ipol, yy);
  }
  return 0;
}
