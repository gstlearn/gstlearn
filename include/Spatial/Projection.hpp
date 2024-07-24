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
#pragma once

#include "gstlearn_export.hpp"

#include "Db/Db.hpp"
#include "Basic/VectorNumT.hpp"

class Db;
class Polygons;

class GSTLEARN_EXPORT Projection
{
public:
  Projection(bool flag_mean, double xcenter=TEST, double ycenter=TEST);
  Projection(bool flag_mean, Db* db);
  Projection(const Projection& r);
  Projection& operator=(const Projection& r);
  virtual ~Projection();

  void operateInPlace(VectorDouble& coor);
  VectorDouble operateInvert(const VectorDouble& coor);
  int operateVecInPlace(VectorDouble& x, VectorDouble& y);
  int operateOnDb(Db *db);
  int operateOnPolygons(Polygons* poly);

  bool isFlagMean() const { return _flagMean; }
  double getXcenter() const { return _xcenter; }
  double getYcenter() const { return _ycenter; }

private:
  bool _flagMean;
  double _xcenter;
  double _ycenter;
};
