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
#include "Basic/AStringable.hpp"

Projection::Projection(bool flag_mean, double xcenter, double ycenter)
    : _flagMean(true),
      _xcenter(xcenter),
      _ycenter(ycenter)
{
}

Projection::Projection(bool flag_mean, Db* db)
: _flagMean(true),
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

  return;
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


int Projection::operateInPlace(VectorDouble& x, VectorDouble& y)
{
  int nech = (int) x.size();
  if (nech != (int) y.size())
  {
    messerr("Arguments 'x' and 'y' should have same dimension");
    return 1;
  }

  // Read the global parameters
  double xc = _xcenter;
  double yc = _ycenter;

  double xx, yy;
  if (_flagMean)
  {
    // Mean projection

    for (int iech = 0; iech < nech; iech++)
    {
      xx = x[iech];
      yy = y[iech];
      x[iech] = (xx - xc) * 60 * cos(yc * GV_PI / 180.);
      y[iech] = (yy - yc) * 60.;
    }
  }
  else
  {
    // Cosine projection

    for (int iech = 0; iech < nech; iech++)
    {
      xx = x[iech];
      yy = y[iech];
      x[iech] = (xx - xc) * 60 * cos((yy - yc + yc * GV_PI) / 180);
      y[iech] = (yy - yc) * 60.;
    }
  }
  return 0;
}

int Projection::db_projection(Db *db)
{
  if (db->getLocNumber(ELoc::X) < 2)
  {
    messerr("This method is dedicated to 2-D space (or more)");
    return 1;
  }

  // Extract the vector of coordinates
  VectorDouble x = db->getCoordinates(0, true);
  VectorDouble y = db->getCoordinates(1, true);

  // Perform the projection
//  if (projectionInPlace(x,  y)) return 1;

  // Store the coordinates (in place)
  db->setCoordinates(0, x, true);
  db->setCoordinates(1, y, true);

  return 0;
}
