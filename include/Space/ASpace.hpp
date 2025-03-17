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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Enum/ESpaceType.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/ICloneable.hpp"
#include <memory>

class SpacePoint;
class Tensor;

/**
  * @brief Base classe for space definitions
  *
  * If the space instance is a sub-space from a parent SpaceComposit
  * _offset is used. Otherwise it is set to 0.
  * Example : if I am RN(1) in RN(2)+RN(1), offset is 2
  */
class GSTLEARN_EXPORT ASpace: public AStringable,
                              public ICloneable
{
protected:
  ASpace(unsigned int ndim);
  ASpace(const ASpace& r);
  ASpace& operator=(const ASpace& r);

public: 
  virtual ~ASpace();

public:
  /// Interface for AStringable
  String toString(const AStringFormat* strfmt = nullptr) const final;

  /// Return the concrete space type
  virtual ESpaceType getType() const = 0;

  ///////////////////////////////////////////////
  // Default behavior that can be overriden

  /// Update the origin of the space
  virtual void setOrigin(const VectorDouble& origin);

  /// Get the number of dimensions
  virtual unsigned int getNDim(int ispace = -1) const;

  /// Get the offset index for coordinates
  virtual unsigned int getOffset(int ispace = -1) const;

  /// Return the space origin coordinates
  virtual const VectorDouble& getOrigin(int ispace = -1) const;

  /// Get the number of space components
  virtual unsigned int getNComponents() const;

  /// Return the space component at index ispace
  virtual std::shared_ptr<const ASpace> getComponent(int ispace = -1) const;

  /// Dump a space in a string (given the space index)
  virtual String toString(const AStringFormat* strfmt, int ispace) const;

  /// Return true if the given space is equal to me (same dimension and space
  /// definition)
  virtual bool isEqual(const ASpace* space) const;

  /// Return all the distances (one by space component) between two space points
  virtual VectorDouble getDistances(const SpacePoint& p1,
                                    const SpacePoint& p2) const;

  virtual void getDistancePointVectInPlace(const SpacePoint& p1,
                                           const std::vector<SpacePoint>& p2,
                                           VectorDouble& res) const
  {
    DECLARE_UNUSED(p1, p2, res)                                      
  };
  ///////////////////////////////////////////////
  /// Not to be overriden
  
  /// Move the given space point by the given vector
  void move(SpacePoint& p1, const VectorDouble& vec) const;

  /// Return the distance between two space points
  double getDistance(const SpacePoint &p1,
                     const SpacePoint &p2,
                     int ispace = -1) const;

  /// Return the distance between two space points with the given tensor
  double getDistance(const SpacePoint& p1,
                     const SpacePoint& p2,
                     const Tensor& tensor,
                     int ispace = -1) const;

  /// Return the distance in frequential domain between two space points with the given tensor
  double getFrequentialDistance(const SpacePoint& p1,
                                const SpacePoint& p2,
                                const Tensor& tensor,
                                int ispace = -1) const;

  /// Return the increment vector between two space points
  VectorDouble getIncrement(const SpacePoint& p1,
                            const SpacePoint& p2,
                            int ispace = -1) const;

  /// Project the coordinates in the given space
  virtual VectorDouble projCoord(const VectorDouble& coord,
                                 int ispace = -1) const;
  
  /// Customize the dimension offset index of the current space
  /// TODO : to be made private
  void setOffset(unsigned int offset) { _offset = offset; }

  static std::shared_ptr<const ASpace> getDefaultSpaceIfNull(const std::shared_ptr<const ASpace>& space);
protected:

  /// Move the given space point by the given vector
  virtual void _move(SpacePoint& p1, const VectorDouble& vec) const = 0;

  /// Return the distance between two space points
  virtual double _getDistance(const SpacePoint& p1,
                              const SpacePoint& p2,
                              int ispace = -1) const = 0;

  /// Return the distance between two space points with the given tensor
  virtual double _getDistance(const SpacePoint& p1,
                              const SpacePoint& p2,
                              const Tensor& tensor,
                              int ispace = -1) const = 0;

  /// Return the distance in frequential domain between two space points with the given tensor
  virtual double _getFrequentialDistance(const SpacePoint& p1,
                                         const SpacePoint& p2,
                                         const Tensor& tensor,
                                         int ispace = -1) const = 0;

  /// Return the increment vector between two space points
  virtual VectorDouble _getIncrement(const SpacePoint& p1,
                                     const SpacePoint& p2,
                                     int ispace = -1) const = 0;
  
  /// Return the increment vector between two space points in a given vector
  virtual void _getIncrementInPlace(const SpacePoint& p1,
                                    const SpacePoint& p2,
                                    VectorDouble& ptemp,
                                    int ispace = -1) const = 0;

protected:
  /// Customize the dimension offset index of the current space
  //void _setOffset(unsigned int offset) { _offset = offset; }
  
protected:
  /// Number of space dimensions
  unsigned int _nDim;
  /// Coordinates of the origin (size = _nDim)
  VectorDouble _origin;
  /// Dimension offset index (0 if single space, relative if sub-space)
  unsigned int _offset;

  /// The next vectors are specified as working members in order to avoid too many allocations
  mutable VectorDouble _work1;
  mutable VectorDouble _work2;

  /// Privilege to SpaceComposit only
  //friend class SpaceComposit; /// TODO : this has no effect (see _setOffset). Why ?
};

typedef std::shared_ptr<const ASpace> ASpaceSharedPtr;
typedef std::vector<ASpaceSharedPtr> ASpaceSharedPtrVector; 