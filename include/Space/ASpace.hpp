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

#include "Enum/ESpaceType.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/ICloneable.hpp"

class SpacePoint;
class Tensor;

class GSTLEARN_EXPORT ASpace : public AStringable, public ICloneable
{
public:
  ASpace(unsigned int ndim);
  ASpace(const ASpace& r);
  ASpace& operator=(const ASpace& r);
  virtual ~ASpace();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Update the origin coordinates
  void setOrigin(const VectorDouble& origin);

  /// Get the number of dimensions
  unsigned int getNDim() const { return _nDim; }

  /// Return the space origin coordinates
  const VectorDouble& getOrigin() const { return _origin; }
  /// Return true if the given space is equal to me
  virtual bool isEqual(const ASpace* space) const;

  /// Return the concrete space type
  virtual ESpaceType getType() const = 0;

  /// Move the given space point by the given vector
  virtual void move(SpacePoint& p1, const VectorDouble& vec) const = 0;

  /// Return the distance between two space points
  virtual double getDistance(const SpacePoint &p1,
                             const SpacePoint &p2) const = 0;

  /// Return the distance between two space points with the given tensor
  virtual double getDistance(const SpacePoint& p1,
                             const SpacePoint& p2,
                             const Tensor& tensor) const = 0;

  /// Return the distance in frequential domain between two space points with the given tensor
  virtual double getFrequentialDistance(const SpacePoint& p1,
                                        const SpacePoint& p2,
                                        const Tensor& tensor) const = 0;

  /// Return the increment vector between two space points
  virtual VectorDouble getIncrement(const SpacePoint& p1,
                                    const SpacePoint& p2) const = 0;

  virtual void _getIncrementInPlace(const SpacePoint &p1,
                                    const SpacePoint &p2,
                                    VectorDouble &ptemp) const { }
  void _getIncrementInPlaceVect(const SpacePoint &p1,
                                const std::vector<SpacePoint> &pv,
                                VectorVectorDouble &res) const;


protected:
  /// Number of space dimensions
  unsigned int _nDim;
  /// Coordinates of the origin (not a space point... ex: sphere center in a long/lat space)
  VectorDouble _origin;

  // The next vectors are specified as working members in order to avoid too many allocations
  mutable VectorDouble _work1;
  mutable VectorDouble _work2;
};

