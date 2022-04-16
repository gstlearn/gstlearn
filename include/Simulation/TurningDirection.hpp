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
#pragma once

#include "gstlearn_export.hpp"
#include "Simulation/ASimulation.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/Vector.hpp"

class Db;

class GSTLEARN_EXPORT TurningDirection
{
public:
  TurningDirection();
  TurningDirection(const TurningDirection& r);
  TurningDirection& operator=(const TurningDirection& r);
  virtual ~TurningDirection();

  const VectorDouble& getAng() const { return _ang; }
  double getAng(int idir) const { return _ang[idir]; }
  double getDXP()   const { return _dxp; }
  double getDYP()   const { return _dyp; }
  double getDZP()   const { return _dzp; }
  double getT00()   const { return _t00; }
  double getTmax()  const { return _tmax; }
  double getTmin()  const { return _tmin; }
  double getScale() const { return _scale; }

  void setAng(const VectorDouble& ang) { _ang = ang; }
  void setAng(int idir, double value) { _ang[idir] = value; }
  void setDXP(double dxp)     { _dxp = dxp; }
  void setDYP(double dyp)     { _dyp = dyp; }
  void setDZP(double dzp)     { _dzp = dzp; }
  void setT00(double t00)     { _t00 = t00; }
  void setTmax(double tmax)   { _tmax = tmax; }
  void setTmin(double tmin)   { _tmin = tmin; }
  void setScale(double scale) { _scale = scale; }

  double projectPoint(const Db* db, int iech) const;
  double projectGrid(const DbGrid* db, int ix, int iy, int iz) const;

private:
  double _tmin; /* Minimum abscissa along line */
  double _tmax; /* Maximum abscissa along line */
  double _scale; /* Scaling factor */
  double _t00; /* Origin along the line */
  double _dxp; /* Increment along X */
  double _dyp; /* Increment along Y */
  double _dzp; /* Increment along Z */
  VectorDouble _ang; /* Angles for the line orientation */
};
