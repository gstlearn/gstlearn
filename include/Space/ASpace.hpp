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

#include "Enum/ESpaceType.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/ICloneable.hpp"

#include <vector>

class SpacePoint;
class Tensor;

class GSTLEARN_EXPORT ASpace : public AStringable, public ICloneable
{
public:
  ASpace(unsigned int ndim, bool addTime = false);
  ASpace(const ASpace& r);
  ASpace& operator=(const ASpace& r);
  virtual ~ASpace();

  /// Return the concrete space type
  virtual ESpaceType getType() const = 0;

  /// Add a space component to me (for exemple RN(1) for time dimension)
  void addSpaceComponent(const ASpace* comp);

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const final;

  /// Update the origin coordinates
  void setOrigin(const VectorDouble& origin);

  /// Get the number of dimensions
  unsigned int getNDim() const;

  /// Get the number of space components
  unsigned int getNComponents() const;

  /// Return the space origin coordinates
  const VectorDouble& getOrigin() const;

  /// Return true if the given space is equal to me (same dimension and space definition)
  bool isEqual(const ASpace* space) const;

  /// Move the given space point by the given vector
  void move(SpacePoint& p1, const VectorDouble& vec) const;

  /// Return the distance between two space points
  double getDistance(const SpacePoint &p1,
                     const SpacePoint &p2) const;

  /// Return the distance between two space points with the given tensor
  double getDistance(const SpacePoint& p1,
                     const SpacePoint& p2,
                     const Tensor& tensor) const;

  /// Return the distance along one direction between two space points
  double getDistance1D(const SpacePoint &p1,
                       const SpacePoint &p2,
                       int idim = 0) const;

  /// Return the distance in frequential domain between two space points with the given tensor
  double getFrequentialDistance(const SpacePoint& p1,
                                const SpacePoint& p2,
                                const Tensor& tensor) const;

  /// Return the increment vector between two space points
  VectorDouble getIncrement(const SpacePoint& p1,
                            const SpacePoint& p2) const;

protected:

  /// Dump a space in a string
  virtual String _toString(const AStringFormat* strfmt, int idx = -1) const;

  /// Return true if the given space is equal to me (same dimension and space definition)
  virtual bool _isEqual(const ASpace* space) const;

  /// Move the given space point by the given vector
  virtual void _move(SpacePoint& p1, const VectorDouble& vec) const = 0;

  /// Return the distance between two space points
  virtual double _getDistance(const SpacePoint &p1,
                              const SpacePoint &p2) const = 0;

  /// Return the distance between two space points with the given tensor
  virtual double _getDistance(const SpacePoint& p1,
                              const SpacePoint& p2,
                              const Tensor& tensor) const = 0;

  /// Return the distance along one direction between two space points
  virtual double _getDistance1D(const SpacePoint &p1,
                                const SpacePoint &p2,
                                int idim = 0) const = 0;

  /// Return the distance in frequential domain between two space points with the given tensor
  virtual double _getFrequentialDistance(const SpacePoint& p1,
                                         const SpacePoint& p2,
                                         const Tensor& tensor) const = 0;

  /// Return the increment vector between two space points
  virtual VectorDouble _getIncrement(const SpacePoint& p1,
                                     const SpacePoint& p2) const = 0;

  /// Return the increment vector between two space points in a given vector
  virtual void _getIncrementInPlace(const SpacePoint &p1, 
                                    const SpacePoint &p2,
                                    VectorDouble &ptemp) const = 0;

protected:
  /// Space composits list
  std::vector<ASpace*> _comps;
  
  /// Number of space dimensions
  unsigned int _nDim;
  /// Coordinates of the origin (not a space point... ex: sphere center in a long/lat space)
  VectorDouble _origin;
  /// Coordinates of the global origin (taking into account composits)
  VectorDouble _globalOrigin;

  // The next vectors are specified as working members in order to avoid too many allocations
  mutable VectorDouble _work1;
  mutable VectorDouble _work2;
};
