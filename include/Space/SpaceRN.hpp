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

#include "Basic/VectorNumT.hpp"
#include "Space/ASpace.hpp"

namespace gstlrn
{
class SpacePoint;
class Tensor;

class GSTLEARN_EXPORT SpaceRN: public ASpace
{
private:
  SpaceRN(size_t ndim);
  SpaceRN(const SpaceRN& r);
  SpaceRN& operator=(const SpaceRN& r);

public:
  virtual ~SpaceRN();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(SpaceRN)

  /// Return the concrete space type
  ESpaceType getType() const override { return ESpaceType::RN; }
  static ASpaceSharedPtr create(Id ndim);
  void getDistancePointVectInPlace(const SpacePoint& p1,
                                   const std::vector<SpacePoint>& p2,
                                   VectorDouble& res,
                                   const VectorInt& ranks) const override;

protected:
  /// Move the given space point by the given vector
  void _move(SpacePoint& p1, const VectorDouble& vec) const override;

  /// Return the distance between two space points
  double _getDistance(const SpacePoint& p1,
                      const SpacePoint& p2,
                      Id ispace = -1) const override;

  /// Return the distance between two space points with the given tensor
  double _getDistance(const SpacePoint& p1,
                      const SpacePoint& p2,
                      const Tensor& tensor,
                      Id ispace = -1) const override;

  /// Return the distance in frequential domain between two space points with the given tensor
  double _getFrequentialDistance(const SpacePoint& p1,
                                 const SpacePoint& p2,
                                 const Tensor& tensor,
                                 Id ispace = -1) const override;

  /// Return the increment vector between two space points for the current space context
  VectorDouble _getIncrement(const SpacePoint& p1,
                             const SpacePoint& p2,
                             Id ispace = -1) const override;

  /// Return the increment vector between two space points in a given vector
  void _getIncrementInPlace(const SpacePoint& p1,
                            const SpacePoint& p2,
                            VectorDouble& ptemp,
                            Id ispace = -1) const override;
};
} // namespace gstlrn