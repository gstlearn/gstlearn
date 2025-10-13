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

#include "Basic/AStringable.hpp"

#include <array>

namespace gstlrn
{
class AShape;
class Db;
class DbGrid;
class ModelBoolean;
class SimuBooleanParam;

class GSTLEARN_EXPORT BooleanObject: public AStringable
{
public:
  BooleanObject(const AShape* shape);
  BooleanObject(const BooleanObject &r);
  BooleanObject& operator=(const BooleanObject &r);
  virtual ~BooleanObject();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  void setMode(Id mode) { _mode = mode; }
  void setOrientation(double orientation) { _orientation = orientation; }
  void setCenter(Id idim, double value) { _center[idim] = value; }
  void setCenter(const VectorDouble& center);
  void setExtension(Id idim, double value) { _extension[idim] = value; }
  void setValue(Id rank, double value) { _values[rank] = value; }

  Id getMode() const { return _mode; }
  double getCenter(Id idim) const { return _center[idim]; }
  double getExtension(Id idim) const { return _extension[idim]; }
  double getOrientation() const { return _orientation; }
  double getValue(Id rank) const { return _values[rank]; }
  const AShape* getToken() const { return _token; }

  static BooleanObject* generate(const DbGrid* dbout,
                          const VectorDouble& cdgrain,
                          const ModelBoolean* tokens,
                          const SimuBooleanParam& boolparam,
                          double eps = EPSILON3);

  bool isCompatiblePore(const Db* db);
  bool isCompatibleGrainAdd(const Db* db);
  bool isCompatibleGrainDelete(const Db* db, Id iptr_cover);
  void projectToGrid(DbGrid* dbout,
                     Id iptr_simu,
                     Id iptr_rank,
                     Id facies,
                     Id rank);
  Id  coverageUpdate(Db* db, Id iptr_cover, Id val);
  VectorDouble getValues() const;

private:
  static bool _invalidTokenFromIntensity(const DbGrid* dbout,
                              const ModelBoolean* tokens,
                              const VectorDouble& coor,
                              double eps = EPSILON3);
  static bool _isPore(const Db* db, Id iech);
  static bool _isGrain(const Db* db, Id iech);
  void _defineBoundingBox(double eps = EPSILON3);
  void _extensionLinkage();
  bool _isInObject(const VectorDouble& coor, Id ndim);

  bool _isInBoundingBox(const VectorDouble& coor, Id ndim);
  static Id  _getCoverageAtSample(const Db* db, Id iptr_cover, Id iech);
  static void _updateCoverageAtSample(Db* db, Id iptr_cover, Id iech, Id ival);
  static void _drawCoordinate(const DbGrid *dbout,
                              const SimuBooleanParam& boolparam,
                              VectorDouble& coor);

private:
  Id _mode;                // 1 for Primary; 2 for Secondary object
  const AShape* _token;     // Token to which the Object belongs
  std::array<double, 3> _center;     // Coordinates of the center of the object
  std::array<double, 3> _extension;  // Extension of the object
  double _orientation;      // Orientation angle for the object (degree)
  std::array<double, 3> _values;     // List of additional arguments
  std::array<std::array<double, 2>, 3> _box;  // Bounding Box containing the object
};
}