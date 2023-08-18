/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
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
class Tensor;

class GSTLEARN_EXPORT SpaceRN : public ASpace {

public:
  SpaceRN(unsigned int ndim);
  SpaceRN(const SpaceRN& r);
  SpaceRN& operator=(const SpaceRN& r);
  virtual ~SpaceRN();

  /// ICloneable interface
  IMPLEMENT_CLONING(SpaceRN)

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static SpaceRN* create(unsigned int ndim);

  /// Return the concrete space type
  ESpaceType getType() const override { return ESpaceType::fromKey("RN"); }
  /// Move the given space point by the given vector
  void move(SpacePoint &p1, const VectorDouble &vec) const override;
  /// Return the distance between two space points
  double getDistance(const SpacePoint &p1, const SpacePoint &p2) const override;
  /// Return the distance between two space points with the given tensor
  double getDistance(const SpacePoint &p1,
                     const SpacePoint &p2,
                     const Tensor &tensor) const override;
  double getFrequentialDistance(const SpacePoint& p1,
                                const SpacePoint& p2,
                                const Tensor& tensor) const override;
  /// Return the increment vector between two space points for the current space context
  VectorDouble getIncrement(const SpacePoint& p1,
                            const SpacePoint& p2) const override;

private:
  void _getIncrementInPlace(const SpacePoint &p1,
                            const SpacePoint &p2,
                            VectorDouble &ptemp) const override;
};
