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
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"

class Db;

class GSTLEARN_EXPORT Projection
{
public:
  Projection(bool flag_mean, double xcenter=TEST, double ycenter=TEST);
  Projection(bool flag_mean, Db* db);
  Projection(const Projection& r);
  Projection& operator=(const Projection& r);
  virtual ~Projection();

  int operateInPlace(VectorDouble& x, VectorDouble& y);

  int db_projection(Db *db);

private:
  bool _flagMean;
  double _xcenter;
  double _ycenter;
};
