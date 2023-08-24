/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "ACalcSimulation.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/VectorNumT.hpp"

class Db;

/**
 * Class for management of Directions used in Turning Band algorithm
 * Remark: The 3-D definition is compulsory (even in 2-D)
 */
class GSTLEARN_EXPORT TurningDirection
{
public:
  TurningDirection();
  TurningDirection(const TurningDirection& r);
  TurningDirection& operator=(const TurningDirection& r);
  virtual ~TurningDirection();

  const VectorDouble& getAng() const { return _ang; }
  double getAng(int i) const { return _ang[i]; }
  double getDXP()   const { return _dxp; }
  double getDYP()   const { return _dyp; }
  double getDZP()   const { return _dzp; }
  double getT00()   const { return _t00; }
  double getTmax()  const { return _tmax; }
  double getTmin()  const { return _tmin; }
  double getScale() const { return _scale; }

  void setAng(const VectorDouble& ang) { _ang = ang; }
  void setAng(int i, double value) { _ang[i] = value; }
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
