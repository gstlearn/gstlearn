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
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"

class ASpace;
class SpacePoint;
class Tensor;

class ASpaceObject : public AStringable
{
public:
  ASpaceObject(const ASpace* space = nullptr);
  ASpaceObject(const ASpaceObject& r);
  ASpaceObject& operator= (const ASpaceObject& r);
  virtual ~ASpaceObject();

  /// Factory for creating the unique global space (optional param can be used for sphere radius)
  static const ASpace* createGlobalSpace(SpaceType type,
                                         unsigned int ndim,
                                         double param = 0.);

  /// Delete the current Global Space (Dangerous function)
  static void destroyGlobalSpace();

  /// Check that the global Space has already been defined
  static bool hasGlobalSpace() { return _globalSpace != nullptr; }

  /// Check the validity of the Space Dimension
  static bool isSpaceDimensionValid(int ndim);

  /// Return the unique global space
  static const ASpace* getGlobalSpace();

  /// Accessor to the current object space context
  const ASpace* getSpace() const { return _space; }
  /// Indicate if I am consistent with my current space context
  bool isConsistent() const { return isConsistent(_space); }
  /// Return unitary vector for the current space context
  VectorDouble getUnitaryVector() const;

  /// Indicate if I am consistent with the provided space
  virtual bool isConsistent(const ASpace* space) const = 0;

  //////////////////////////////////////////////////////////
  /// Shortcuts to ASpace methods

  /// Return the number of dimension of the current space context
  unsigned int getNDim() const;
  /// Return the current space context origin
  const SpacePoint& getOrigin() const;
  /// Return the distance between two space points for the current space context
  double getDistance(const SpacePoint& p1, const SpacePoint& p2) const;
  /// Return the increment vector between two space points for the current space context
  VectorDouble getIncrement(const SpacePoint& p1, const SpacePoint& p2) const;

private:
  /// Unique and global space
  static ASpace* _globalSpace;
  /// True if Global Space has been created by default
  static bool _fakeSpace;
  /// Current space context
  const ASpace* _space;
};
