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

#include "Space/ASpace.hpp"
#include "Basic/VectorNumT.hpp"

class SpacePoint;

class GSTLEARN_EXPORT SpaceSN: public ASpace
{

public:
  SpaceSN(unsigned int ndim, double radius);
  SpaceSN(const SpaceSN &r);
  SpaceSN& operator=(const SpaceSN &r);
  virtual ~SpaceSN();

  /// ICloneable interface
IMPLEMENT_CLONING(SpaceSN)

  virtual String toString(const AStringFormat *strfmt = nullptr) const override;

  /// Return true if the given space is equal to me
  virtual bool isEqual(const ASpace *space) const override;

  /// Return the concrete space type
  ESpaceType getType() const override
  {
    return ESpaceType::fromKey("SN");
  }

  /// Move the given space point by the given vector
  void move(SpacePoint &p1, const VectorDouble &vec) const override;
  /// Return the distance between two space points
  double getDistance(const SpacePoint &p1, const SpacePoint &p2) const override;
  /// Return the distance between two space points with the given tensor
  double getDistance(const SpacePoint &p1,
                     const SpacePoint &p2,
                     const Tensor &tensor) const override;
  /// Return the distance along one direction between two space points
  double getDistance1D(const SpacePoint &p1,
                       const SpacePoint &p2,
                       int idim = 0) const override;

  /// Return the distance in frequential domain between two space points with the given tensor
  double getFrequentialDistance(const SpacePoint &p1,
                                const SpacePoint &p2,
                                const Tensor &tensor) const override;

  /// Return the increment vector between two space points for the current space context
  VectorDouble getIncrement(const SpacePoint &p1, const SpacePoint &p2) const
      override;

  double getRadius() const
  {
    return _radius;
  }

private:
  /// Sphere radius
  double _radius;
};

