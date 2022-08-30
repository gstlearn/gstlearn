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
#include "Space/ASpace.hpp"
#include "Basic/Vector.hpp"

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

  /// Return the concrete space type
  ESpaceType getType() const override { return ESpaceType::SPACE_RN; }
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
  double getFrequentialDistance(const SpacePoint& p1,
                                const SpacePoint& p2,
                                const Tensor& tensor) const override;
  /// Return the increment vector between two space points for the current space context
  VectorDouble getIncrement(const SpacePoint& p1,
                            const SpacePoint& p2) const override;
};
