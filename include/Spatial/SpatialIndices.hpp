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
#include "geoslib_define.h"

#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"
#include "Spatial/Projection.hpp"

class Db;
class MatrixRectangular;

class GSTLEARN_EXPORT SpatialIndices : public AStringable
{
public:
  SpatialIndices(Db* db);
  SpatialIndices(const SpatialIndices &r);
  SpatialIndices &operator=(const SpatialIndices &r);
  virtual ~SpatialIndices();

  /// Interface to AStringable
  virtual String toString(const AStringFormat *strfmt = nullptr) const override;

  void cgiPrint() const;
  int computeCGI(Db *db, const String &name);
  void spatial(Db *db);

  MatrixRectangular getMatrixEllipse() const;
  MatrixRectangular getMatrixInertia() const;
  VectorVectorDouble getAxes() const;

private:
  int _getData(bool flag_w, int iech, const String &name, VectorDouble &coor,
               double *wz);

private:
  Db *_db; // Not to be deleted
  VectorDouble _center; // Vector for the center of gravity
  VectorDouble _mvalues;       // Vector of eigen values (normalized)
  MatrixRectangular _mvectors; // Array of eigen vectors
  double _inertia; // Value of the inertia
  double _wztot;   // Sum of weights
  double _iso;     // Iso index
  int    _nvalid;     // Number of valid samples

  double _theta;
  double _ra;
  double _rb;

  double _totab;
  double _parea;
  double _eqarea;
};

