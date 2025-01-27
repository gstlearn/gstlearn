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
#include "Enum/ESpaceType.hpp"

#include <memory>
#include <vector>

class SpacePoint;
class Tensor;

class GSTLEARN_EXPORT SpaceComposite : public ASpace
{
private:
  SpaceComposite(const std::vector<ASpaceSharedPtr>& vectspace = std::vector<ASpaceSharedPtr>());
  SpaceComposite(const SpaceComposite& r);
  SpaceComposite& operator=(const SpaceComposite& r);
  
public:
  virtual ~SpaceComposite();

  /// ICloneable interface
  IMPLEMENT_CLONING(SpaceComposite)

  static std::shared_ptr<SpaceComposite> create(const ASpaceSharedPtrVector& vectspace = ASpaceSharedPtrVector());
  /// Return the concrete space type
  ESpaceType getType() const override { return ESpaceType::COMPOSITE; }

  /// Update the origin of the space
  void setOrigin(const VectorDouble& origin) override;

  /// Get the number of dimensions
  unsigned int getNDim(int ispace = -1) const override;

  /// Get the offset index for coordinates
  unsigned int getOffset(int ispace = -1) const override;
  
  /// Return the space origin coordinates
  const VectorDouble& getOrigin(int ispace = -1) const override;

  /// Get the number of space components
  unsigned int getNComponents() const override;

  /// Return the space component at index ispace
  ASpaceSharedPtr getComponent(int ispace = -1) const override;

  /// Dump a space in a string (given the space index)
  String toString(const AStringFormat* strfmt, int ispace) const override;

  /// Return true if the given space is equal to me (same dimension and space
  /// definition)
  bool isEqual(const ASpace* space) const override;

  /// Return all the distances (one by space component) between two space points
  VectorDouble getDistances(const SpacePoint& p1,
                            const SpacePoint& p2) const override;
  
  /////////////////////////////////////////////
  
  /// Add a space component to me (for exemple RN(1) for time dimension)
  void addSpaceComponent(const ASpaceSharedPtr& comp);

protected:

  /// Move the given space point by the given vector
  void _move(SpacePoint& p1, const VectorDouble& vec) const override;

  /// Return the distance between two space points
  double _getDistance(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ispace = -1) const override;

  /// Return the distance between two space points with the given tensor
  double _getDistance(const SpacePoint& p1,
                      const SpacePoint& p2,
                      const Tensor& tensor,
                      int ispace = -1) const override;

  /// Return the distance in frequential domain between two space points with
  /// the given tensor
  double _getFrequentialDistance(const SpacePoint& p1,
                                 const SpacePoint& p2,
                                 const Tensor& tensor,
                                 int ispace = -1) const override;

  /// Return the increment vector between two space points for the current space
  /// context
  VectorDouble _getIncrement(const SpacePoint& p1,
                             const SpacePoint& p2,
                             int ispace = -1) const override;

  /// Return the increment vector between two space points in a given vector
  void _getIncrementInPlace(const SpacePoint& p1,
                            const SpacePoint& p2,
                            VectorDouble& ptemp,
                            int ispace = -1) const override;


private:
  /// Space composits list
  std::vector<std::shared_ptr<ASpace> > _comps;
};
