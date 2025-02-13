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

class GSTLEARN_EXPORT SpacePoint : public ASpaceObject
{
public:
  SpacePoint(const ASpaceSharedPtr& space = ASpaceSharedPtr());
  SpacePoint(const SpacePoint& r);
  SpacePoint(const VectorDouble& coord, int iech = -1,
             const ASpaceSharedPtr& space = ASpaceSharedPtr());
  SpacePoint& operator=(const SpacePoint& r);
  virtual ~SpacePoint();

  SpacePoint spacePointOnSubspace(int ispace = -1) const;

  /// TODO : should aslo test the space definition
  bool operator==(const SpacePoint& v) const { return (_coord == v._coord); }

  constvect getCoords() const;
  
  vect getCoordRef() { return vect(_coord.data(), getNDim()); }
  double getCoord(int idim) const; 
  void setCoord(double coord);
  void setCoord(int i, double val) { _coord[i] = val; }
  void setCoords(const VectorDouble& coord);
  void setCoords(const double* coord, int size);
  void setIech(int iech) const { _iech = iech; }
  int getIech() const { return _iech; }
  void setTarget(bool target) const { _target = target; }
  bool isTarget() const { return _target; }
  /// Return true if the point is consistent with the provided space
  virtual bool isConsistent(const ASpace* space) const override;

  /// Move me by the given vector
  void move(const VectorDouble& vec);
  /// Return the distance between me and another point
  double getDistance(const SpacePoint& pt, int ispace = -1) const;
  /// Return all the distance (space composits) between me and another point
  VectorDouble getDistances(const SpacePoint& pt) const;
  /// Return the increment vector between me and another point
  VectorDouble getIncrement(const SpacePoint& pt, int ispace = -1) const;
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
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

protected:
  /// Points coordinates (whatever the space context)
  VectorDouble _coord;
  mutable int  _iech;
  mutable bool _target;
};
