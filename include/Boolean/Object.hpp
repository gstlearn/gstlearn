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
#include "geoslib_define.h"

#include "Basic/AStringable.hpp"

class AToken;
class Db;
class DbGrid;
class Tokens;

class GSTLEARN_EXPORT Object: public AStringable
{
public:
  Object(const AToken* atoken);
  Object(const Object &r);
  Object& operator=(const Object &r);
  virtual ~Object();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void setMode(int mode) { _mode = mode; }
  void setOrientation(double orientation) { _orientation = orientation; }
  void setCenter(int idim, double value) { _center[idim] = value; }
  void setCenter(const VectorDouble& center) { _center = center; }
  void setExtension(int idim, double value) { _extension[idim] = value; }
  void setValue(int rank, double value) { _values[rank] = value; }

  int getMode() const { return _mode; }
  double getCenter(int idim) const { return _center[idim]; }
  double getExtension(int idim) const { return _extension[idim]; }
  double getOrientation() const { return _orientation; }
  double getValue(int rank) const { return _values[rank]; }
  const AToken* getToken() const { return _token; }

  static Object* generate(const DbGrid* dbout,
                          const VectorDouble& cdgrain,
                          const Tokens* tokens,
                          bool flagStat,
                          double thetaCst,
                          const VectorDouble& dilate,
                          int maxiter = 100000,
                          double eps = EPSILON3);

  void blank();
  bool isCompatiblePore(const Db* db);
  bool isCompatibleGrainAdd(const Db* db);
  bool isCompatibleGrainDelete(const Db* db);
  void projectToGrid(DbGrid* dbout,
                     int iptr_simu,
                     int iptr_rank,
                     int facies,
                     int rank);
  int  coverageUpdate(Db* db, int val);

private:
  static bool _checkIntensity(const DbGrid* dbout,
                              const VectorDouble& coor,
                              bool flagStat,
                              double thetaCst,
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
                              const VectorDouble& dilate,
                              VectorDouble& coor);

private:
  int _mode;                // 1 for Primary; 2 for Secondary object
  const AToken* _token;     // Token to which the Object belongs
  VectorDouble _center;     // Coordinates of the center of the object
  VectorDouble _extension;  // Extension of the object
  double _orientation;      // Orientation angle for the object (degree)
  VectorDouble _values;     // List of additional arguments
  VectorVectorDouble _box;  // Bounding Box containing the object
};
