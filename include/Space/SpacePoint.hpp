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

#include "Space/ASpace.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Space/ASpaceObject.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT SpacePoint : public ASpaceObject
{
public:
  SpacePoint(const ASpaceSharedPtr& space = ASpaceSharedPtr());
  SpacePoint(const SpacePoint& r);
  SpacePoint(const VectorDouble& coord, Id iech = -1,
             const ASpaceSharedPtr& space = ASpaceSharedPtr());
  SpacePoint& operator=(const SpacePoint& r);
  virtual ~SpacePoint();

  SpacePoint spacePointOnSubspace(Id ispace = -1) const;

  /// TODO : should also test the space definition
  bool operator==(const SpacePoint& v) const { return (_coord == v._coord); }

#ifndef SWIG
  constvect getCoordsView() const { return {_coord.data(), getNDim()}; }
  vect getCoordsView() { return {_coord.data(), getNDim()}; }
#endif

  const VectorDouble& getCoords() const { return _coord; }
  VectorDouble& getCoordsUnprotected() { return _coord; }

  double getCoord(Id idim) const; 
  void setCoord(double coord);
  void setCoord(Id i, double val) { _coord[i] = val; }
  void setCoords(const VectorDouble& coord);

  void setCoords(const double* coord, Id size);
  void setIech(Id iech) const { _iech = iech; }
  void setProjected(bool status) { _isProjected = status; }
  Id getIech() const { return _iech; }
  bool isProjected() const { return _isProjected; }
  /// Return true if the point is consistent with the provided space
  bool isConsistent(const ASpace* space) const override;

  /// Move me by the given vector
  void move(const VectorDouble& vec);
  /// Return the distance between 'this' and another point
  double getDistance(const SpacePoint& pt, Id ispace = -1) const;
  /// Return all the distance (space composits) between 'this' and another point
  VectorDouble getDistances(const SpacePoint& pt) const;
  /// Return the increment vector between 'this' and another point
  VectorDouble getIncrement(const SpacePoint& pt, Id ispace = -1) const;
  void getIncrementInPlace(VectorDouble &inc, const SpacePoint& pt, Id ispace = -1) const;
  /// Fill with TEST values to simulate a missing Space Point
  void setFFFF();
  /// Check if the SpacePoint is actually defined
  bool isFFFF() const;
  /// Return the cosine of the angle between the bipoint and a reference direction
  double getCosineToDirection(const SpacePoint &T2,
                              const VectorDouble &codir) const;
  /// Return the orthogonal distance between a bipair and a reference direction
  double getOrthogonalDistance(const SpacePoint &P2,
                               const VectorDouble &codir) const;

  /// Initialize coordinates from angles /// TODO : to be removed
  void setCoordFromAngle(const VectorDouble& angles);

  /// Convert space point to string
  String toString(const AStringFormat* strfmt = nullptr) const override;

protected:
  /// Points coordinates (whatever the space context)
  VectorDouble _coord; // Coordinates (initial or projected)
  mutable VectorDouble _delta; // Increment results
  mutable Id _iech;   // Absolute rank of the sample within the Db
  mutable bool _isProjected; // True if the coordinates are projected
};
}