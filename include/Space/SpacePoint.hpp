/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Space/ASpaceObject.hpp"

//#include "Basic/VectorT.hpp"

class GSTLEARN_EXPORT SpacePoint : public ASpaceObject
{
public:
  SpacePoint(const ASpace* space = nullptr);
  SpacePoint(const SpacePoint& r);
  SpacePoint(const VectorDouble& coord,
             const ASpace* space = nullptr);
  SpacePoint& operator=(const SpacePoint& r);
  virtual ~SpacePoint();

  bool operator==(const SpacePoint& v) const { return (_coord == v._coord); }

  const VectorDouble& getCoord() const { return _coord; }
  double getCoord(int idim) const { return _coord[idim]; }

  void setCoord(double coord);
  void setCoord(const VectorDouble& coord);
  void setCoord(int i, double val){_coord[i] = val;}

  /// Return true if the point is consistent with the provided space
  virtual bool isConsistent(const ASpace* space) const override;

  /// Move me by the given vector
  void move(const VectorDouble& vec);
  /// Return the distance between me and another point
  double getDistance(const SpacePoint& pt) const;
  /// Return the increment vector between me and another point
  VectorDouble getIncrement(const SpacePoint& pt) const;

  /// Initialize coordinates from angles /// TODO : to be removed
  void setCoordFromAngle(const VectorDouble& angles);

  /// Convert space point to string
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

protected:
  /// Points coordinates (whatever the space context)
  VectorDouble _coord;
};
