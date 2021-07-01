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

#include "Space/ASpace.hpp"

#include "Basic/Vector.hpp"

class SpacePoint;

class SpaceSN : public ASpace {

public:
  SpaceSN(unsigned int ndim, double radius);
  SpaceSN(const SpaceSN& r);
  SpaceSN& operator=(const SpaceSN& r);
  virtual ~SpaceSN();

  /// Return true if the given space is equal to me
  virtual bool isEqual(const ASpace* space) const override;

  /// Return the concrete space type
  SpaceType getType() const override { return SPACE_SN; }

  /// Move the given space point by the given vector
  void move(SpacePoint& p1,
            const VectorDouble& vec) const override;
  /// Return the distance between two space points
  double getDistance(const SpacePoint& p1,
                     const SpacePoint& p2) const override;
  /// Return the distance between two space points with the given tensor
  double getDistance(const SpacePoint& p1,
                     const SpacePoint& p2,
                     const Tensor& tensor) const override;

  /// Return the increment vector between two space points for the current space context
  VectorDouble getIncrement(const SpacePoint& p1,
                            const SpacePoint& p2) const override;

private:
  /// Sphere radius
  double _radius;
};

