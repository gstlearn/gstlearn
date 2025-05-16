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
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void setMode(int mode) { _mode = mode; }
  void setOrientation(double orientation) { _orientation = orientation; }
  void setCenter(int idim, double value) { _center[idim] = value; }
  void setCenter(const VectorDouble& center);
  void setExtension(int idim, double value) { _extension[idim] = value; }
  void setValue(int rank, double value) { _values[rank] = value; }

  int getMode() const { return _mode; }
  double getCenter(int idim) const { return _center[idim]; }
  double getExtension(int idim) const { return _extension[idim]; }
  double getOrientation() const { return _orientation; }
  double getValue(int rank) const { return _values[rank]; }
  const AShape* getToken() const { return _token; }

  static BooleanObject* generate(const DbGrid* dbout,
                          const VectorDouble& cdgrain,
                          const ModelBoolean* tokens,
                          const SimuBooleanParam& boolparam,
                          double eps = EPSILON3);

  bool isCompatiblePore(const Db* db);
  bool isCompatibleGrainAdd(const Db* db);
  bool isCompatibleGrainDelete(const Db* db, int iptr_cover);
  void projectToGrid(DbGrid* dbout,
                     int iptr_simu,
                     int iptr_rank,
                     int facies,
                     int rank);
  int  coverageUpdate(Db* db, int iptr_cover, int val);
  VectorDouble getValues() const;

private:
  static bool _invalidTokenFromIntensity(const DbGrid* dbout,
                              const ModelBoolean* tokens,
                              const VectorDouble& coor,
                              double eps = EPSILON3);
  static bool _isPore(const Db* db, int iech);
  static bool _isGrain(const Db* db, int iech);
  void _defineBoundingBox(double eps = EPSILON3);
  void _extensionLinkage();
  bool _isInObject(const VectorDouble& coor, int ndim);

  bool _isInBoundingBox(const VectorDouble& coor, int ndim);
  static int  _getCoverageAtSample(const Db* db, int iptr_cover, int iech);
  static void _updateCoverageAtSample(Db* db, int iptr_cover, int iech, int ival);
  static void _drawCoordinate(const DbGrid *dbout,
                              const SimuBooleanParam& boolparam,
                              VectorDouble& coor);

private:
  int _mode;                // 1 for Primary; 2 for Secondary object
  const AShape* _token;     // Token to which the Object belongs
  std::array<double, 3> _center;     // Coordinates of the center of the object
  std::array<double, 3> _extension;  // Extension of the object
  double _orientation;      // Orientation angle for the object (degree)
  std::array<double, 3> _values;     // List of additional arguments
  std::array<std::array<double, 2>, 3> _box;  // Bounding Box containing the object
};
