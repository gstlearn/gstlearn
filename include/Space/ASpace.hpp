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
  /// Note: The given argument is cloned
  void addSpaceComponent(const ASpace* comp);

  /// Get the number of space components
  unsigned int getNComponents() const;

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const final;

  /// Update the origin of the space (must take into account composits)
  void setOrigin(const VectorDouble& origin);

  /// Get the number of dimensions (if ispace is negative, return the global number of )
  unsigned int getNDim(int ispace = -1) const;

  /// Return the space origin coordinates
  const VectorDouble& getOrigin(int ispace = -1) const;

  /// Return the dimension offset index
  unsigned int getDimOffset() const { return _iDimOffset; }

  /// Return true if the given space is equal to me (same dimension and space definition)
  bool isEqual(const ASpace* space) const;

  /// Move the given space point by the given vector
  void move(SpacePoint& p1, const VectorDouble& vec) const;

  /// Return the distance between two space points
  double getDistance(const SpacePoint &p1,
                     const SpacePoint &p2,
                     int ispace = 0) const;

  /// Return the distance between two space points with the given tensor
  double getDistance(const SpacePoint& p1,
                     const SpacePoint& p2,
                     const Tensor& tensor,
                     int ispace = 0) const;

  /// Return all the distances (one by space component) between two space points
  VectorDouble getDistances(const SpacePoint &p1,
                            const SpacePoint &p2) const;

  /// Return the distance along one direction between two space points
  double getDistance1D(const SpacePoint &p1,
                       const SpacePoint &p2,
                       unsigned int idim = 0) const;

  /// Return the distance in frequential domain between two space points with the given tensor
  double getFrequentialDistance(const SpacePoint& p1,
                                const SpacePoint& p2,
                                const Tensor& tensor,
                                int ispace = 0) const;

  /// Return the increment vector between two space points
  VectorDouble getIncrement(const SpacePoint& p1,
                            const SpacePoint& p2,
                            int ispace = 0) const;

protected:

  /// Internal usage only
  void _setDimOffset(unsigned int idim) { _iDimOffset = idim; }

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

  /// Return the distance for a given pair of coordinates along one direction
  virtual double _getDistance1D(double c1, double c2) const = 0;

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
  /// Number of space dimensions (not taking into account composits)
  unsigned int _nDim;
  /// Coordinates of the origin (not taking into account composits)
  VectorDouble _origin;
  /// Dimension offset index (for space composit)
  unsigned int _iDimOffset;

  /// Space composits list
  std::vector<ASpace *> _comps;
  /// Numnber of space dimensions  (taking into account composits)
  unsigned int _globalNDim;
  /// Coordinates of the global origin (taking into account composits)
  VectorDouble _globalOrigin;

  // The next vectors are specified as working members in order to avoid too many allocations
  mutable VectorDouble _work1;
  mutable VectorDouble _work2;
};
