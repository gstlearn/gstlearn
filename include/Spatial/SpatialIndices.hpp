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

  int  computeCGI(const String &name);
  void spatial(const String &name);
  double getLIC(const String &name1, const String &name2);
  double getGIC(const String &name1, const String &name2);

  MatrixRectangular getMatrixEllipse() const;
  MatrixRectangular getMatrixInertia() const;
  VectorVectorDouble getAxes() const;
  VectorDouble getAxe(int rank) const;
  const VectorDouble& getCenter() const { return _center; }
  double getInertia() const { return _inertia; }
  double getIso() const { return _iso; }
  VectorVectorDouble getQT(const String &name) const;
  double getMicroStructure(const String &name, double h0,
                           const Polygons *polygon = nullptr,
                           double dlim = 0., int ndisc = 100);
  std::vector<SpacePoint> getPatches(const String &name, double Dmin, double Amin = 0);

private:
  bool _discardData(bool flag_w, int iech, const String &name,
                    VectorDouble &coor, double *value, double *weight,
                    double *wz) const;

private:
  Db *_db;                     // Not to be deleted
  VectorDouble _center;        // Vector for the center of gravity
  VectorDouble _mvalues;       // Vector of eigen values (normalized)
  MatrixRectangular _mvectors; // Array of eigen vectors
  double _inertia;             // Value of the inertia
  double _wztot;               // Sum of weights
  double _iso;                 // Iso index
  int _nvalid;                 // Number of valid samples

  double _theta;
  double _ra;
  double _rb;
  };
