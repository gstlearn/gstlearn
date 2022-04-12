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
#include "Space/Space.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"

class ASpace;
class SpacePoint;

/**
 * This class is the base class for all objects that need to know what is its space definition.
 * All ASpaceObject can access to the number of space dimensions and can ask to calculate
 * a distance between two ASpaceObjects.
 *
 * This class also stores a unique (static) default global space that will be used as default space
 * when creating a ASpaceObject (without a predefined space). It is possible to modify the default
 * space definition at any time. Space definition of pre-existing ASpaceObjects remains the same.
 * (no more shared pointer)
 */
class GSTLEARN_EXPORT ASpaceObject : public AStringable
{
public:
  ASpaceObject(const ASpace* space = nullptr);
  ASpaceObject(const ASpace& space);
  ASpaceObject(const ASpaceObject& r);
  ASpaceObject& operator= (const ASpaceObject& r);
  virtual ~ASpaceObject();

  /// AStringable interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// (Re)Defining the unique default global space
  static void defineDefaultSpace(SpaceType type,
                                 unsigned int ndim,
                                 double param = 0.);
  /// Return a clone of the unique default global space
  static const ASpace* cloneDefaultSpace();

public:
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

  /// Return the current space context origin coordinates
  const VectorDouble& getOrigin() const;

  /// Return the distance between two space points for the current space context
  double getDistance(const SpacePoint& p1, const SpacePoint& p2) const;

  /// Return the increment vector between two space points for the current space context
  VectorDouble getIncrement(const SpacePoint& p1, const SpacePoint& p2) const;

protected:
  /// Current space context of the object
  const ASpace* _space;

private:
  /// Unique default global space
  static ASpace* _defaultSpace;
};
