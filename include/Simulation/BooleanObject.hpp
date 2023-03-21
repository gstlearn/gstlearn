/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Basic/AStringable.hpp"

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
  bool isCompatibleGrainDelete(const Db* db);
  void projectToGrid(DbGrid* dbout,
                     int iptr_simu,
                     int iptr_rank,
                     int facies,
                     int rank);
  int  coverageUpdate(Db* db, int val);
  VectorDouble getValues() const;

private:
  static bool _checkIntensity(const DbGrid* dbout,
                              const ModelBoolean* tokens,
                              const VectorDouble& coor,
                              double eps = EPSILON3);
  bool _isPore(const Db* db, int iech);
  bool _isGrain(const Db* db, int iech);
  void _defineBoundingBox(double eps = EPSILON3);
  void _extensionLinkage();
  bool _checkObject(const VectorDouble& coor, int ndim);

  bool _checkBoundingBox(const VectorDouble& coor, int ndim);
  int  _getCoverageAtSample(const Db* db, int iech);
  void _updateCoverageAtSample(Db* db, int iech, int ival);
  static void _drawCoordinate(const DbGrid *dbout,
                              const SimuBooleanParam& boolparam,
                              VectorDouble& coor);

private:
  int _mode;                // 1 for Primary; 2 for Secondary object
  const AShape* _token;     // Token to which the Object belongs
  VectorDouble _center;     // Coordinates of the center of the object
  VectorDouble _extension;  // Extension of the object
  double _orientation;      // Orientation angle for the object (degree)
  VectorDouble _values;     // List of additional arguments
  VectorVectorDouble _box;  // Bounding Box containing the object
};
