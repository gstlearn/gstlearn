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

#include "Basic/Vector.hpp"

#include "Space/ASpaceObject.hpp"

//#include "Basic/VectorT.hpp"

class SpacePoint : public ASpaceObject
{
public:
  SpacePoint(const ASpace* space = nullptr);
  SpacePoint(const VectorDouble& coord,
             const ASpace* space = nullptr);
  virtual ~SpacePoint();

  bool operator==(const SpacePoint& v) const { return (_coord == v._coord); }

  const VectorDouble& getCoord() const { return _coord; }

  void setCoord(double coord);
  void setCoord(const VectorDouble& coord);

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
  virtual std::string toString(int level = 0) const override;

protected:
  /// Points coordinates (whatever the space context)
  VectorDouble _coord;
};
