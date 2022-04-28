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
  String toString(const AStringFormat* /*strfmt*/) const;

  void setOrientation(double orientation) { _orientation = orientation; }
  void setCenter(int idim, double value) { _center[idim] = value; }
  void setExtension(int idim, double value) { _extension[idim] = value; }
  void setValue(int rank, double value) { _values[rank] = value; }

  double getCenter(int idim) const { return _center[idim]; }
  double getExtension(int idim) const { return _extension[idim]; }
  double getOrientation() const { return _orientation; }
  double getValue(int rank) const { return _values[rank]; }

  void blank();
  Object* generateObject(const DbGrid* dbout,
                         const VectorDouble& cdgrain,
                         Tokens* tokens,
                         int maxiter,
                         bool flagStat,
                         double thetaCst,
                         double eps = EPSILON3);
  void projectToGrid(DbGrid* dbout,
                     int iptr_simu,
                     int iptr_rank,
                     int facies,
                     int rank);

private:
  bool _isPore(const Db* db, int iech);
  bool _isGrain(const Db* db, int iech);
  void _defineBoundingBox(double eps = EPSILON3);
  void _extensionLinkage();
  bool _checkIntensity(const DbGrid* dbout,
                       const VectorDouble& coor,
                       bool flagStat,
                       double thetaCst,
                       double eps = EPSILON3);
  bool _checkObject(const VectorDouble& coor, int ndim);
  bool _checkPore(const Db* db);
  bool _checkBoundingBox(const VectorDouble& coor, int ndim);
  bool _canBeAdded(const Db* db);
  bool _canBeDeleted(const Db* db);
  int  _getCoverageAtSample(const Db* db, int iech);
  void _updateCoverageAtSample(Db* db, int iech, int ival);
  int  _coverageUpdate(Db* db, int val);

private:
  const AToken* _token;     // Token to which the Object belongs
  VectorDouble _center;     // Coordinates of the center of the object
  VectorDouble _extension;  // Extension of the object
  double _orientation;      // Orientation angle for the object (degree)
  VectorDouble _values;     // List of additional arguments
  VectorVectorDouble _box;  // Bounding Box containing the object
};
