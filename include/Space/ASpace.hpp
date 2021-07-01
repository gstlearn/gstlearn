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

#include "Space/Space.hpp"
#include "Space/SpacePoint.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"

class SpacePoint;
class Tensor;

class ASpace : public AStringable
{
public:
  ASpace(unsigned int ndim);
  ASpace(const ASpace& r);
  ASpace& operator=(const ASpace& r);
  virtual ~ASpace();

  virtual std::string toString(int level = 0) const override;

  /// Update the origin
  void setOrigin(const SpacePoint& origin);

  /// Get the number of dimensions
  unsigned int getNDim() const { return _nDim; }
  /// Return the space origin
  const SpacePoint& getOrigin() const { return _origin; }

  /// Return true if the given space is equal to me
  virtual bool isEqual(const ASpace* space) const;

  /// Return the concrete space type
  virtual SpaceType getType() const = 0;

  /// Move the given space point by the given vector
  virtual void move(SpacePoint& p1,
                    const VectorDouble& vec) const = 0;
  /// Return the distance between two space points
  virtual double getDistance(const SpacePoint& p1,
                             const SpacePoint& p2) const = 0;

  /// Return the distance between two space points with the given tensor
  virtual double getDistance(const SpacePoint& p1,
                             const SpacePoint& p2,
                             const Tensor& tensor) const = 0;

  /// Return the increment vector between two space points
  virtual VectorDouble getIncrement(const SpacePoint& p1,
                                    const SpacePoint& p2) const = 0;

protected:
  /// Number of space dimensions
  unsigned int _nDim;
  /// Origin of the space
  SpacePoint _origin;
};

