/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "Simulation/TurningDirection.hpp"
#include "Db/DbGrid.hpp"

TurningDirection::TurningDirection()
    : _tmin(0.),
      _tmax(0.),
      _scale(0.),
      _t00(0.),
      _dxp(0.),
      _dyp(0.),
      _dzp(0.),
      _ang()
{
  _ang.resize(3);
}

TurningDirection::TurningDirection(const TurningDirection &r)
    : _tmin(r._tmin),
      _tmax(r._tmax),
      _scale(r._scale),
      _t00(r._t00),
      _dxp(r._dxp),
      _dyp(r._dyp),
      _dzp(r._dzp),
      _ang(r._ang)
{
}

TurningDirection& TurningDirection::operator=(const TurningDirection &r)
{
  if (this != &r)
  {
    _tmin = r._tmin;
    _tmax = r._tmax;
    _scale = r._scale;
    _t00 = r._t00;
    _dxp = r._dxp;
    _dyp = r._dyp;
    _dzp = r._dzp;
    _ang = r._ang;
  }
  return *this;
}

TurningDirection::~TurningDirection()
{
}

/*****************************************************************************/
/*!
 **  Calculates the projection of a grid node on a turning band
 **
 ** \return  Projection value
 **
 ** \param[in]  db      Db structure
 ** \param[in]  ix      grid index along X
 ** \param[in]  iy      grid index along Y
 ** \param[in]  iz      grid index along Z
 **
 *****************************************************************************/
double TurningDirection::projectGrid(const DbGrid* db,
                                     int ix,
                                     int iy,
                                     int iz) const
{
  double xyz[3];
  int indg[3];

  indg[0] = ix;
  indg[1] = iy;
  indg[2] = iz;
  grid_to_point(db, indg, nullptr, xyz);

  double t = 0.;
  for (int idim = 0; idim < db->getNDim(); idim++)
    t += xyz[idim] * _ang[idim];
  return t;
}

/*****************************************************************************/
/*!
 **  Calculates the projection of a point on a turning band
 **
 ** \return  Projection value
 **
 ** \param[in]  db      Db structure
 ** \param[in]  iech    rank of the sample
 **
 *****************************************************************************/
double TurningDirection::projectPoint(const Db *db, int iech) const
{
  double t = 0.;
  for (int idim = 0; idim < db->getNDim(); idim++)
    t += db->getCoordinate(iech, idim) * _ang[idim];
  return (t);
}
