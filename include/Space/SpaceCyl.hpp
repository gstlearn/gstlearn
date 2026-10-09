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

  /**
   * @brief Cylindrical surface embedded in R3.
   *
   * A point is stored using its 3 cartesian coordinates (x,y,z).
   * The cylinder is defined by:
   *  - a point through which the axis passes (_axisOrigin), default (0,0,0)
   *  - a unit direction vector for the axis (_axisDir), default (1,0,0)
   *  - a radius (_radius)
   *
   * All points are assumed to lie exactly on the cylinder surface
   * (i.e. the norm of their radial component is equal to _radius).
   */
  class GSTLEARN_EXPORT SpaceCyl: public ASpace
  {
  private:
    SpaceCyl(
      size_t ndim,
      const VectorDouble& axisOrigin,
      const VectorDouble& axisDir,
      double radius);
    SpaceCyl(const SpaceCyl& r);
    SpaceCyl& operator=(const SpaceCyl& r);

  public:
    virtual ~SpaceCyl();

    /// ICloneable interface
    IMPLEMENT_CLONING(SpaceCyl)

    /// Return the concrete space type
    ESpaceType getType() const override { return ESpaceType::CYL; }

    static ASpaceSharedPtr create(
      Id ndim = 3,
      const VectorDouble& axisOrigin = {0., 0., 0.},
      const VectorDouble& axisDir = {1., 0., 0.},
      double radius = 1.);

    /// Return the cylinder radius
    double getRadius() const { return _radius; }

    /// Return the point through which the axis passes
    const VectorDouble& getAxisOrigin() const { return _axisOrigin; }

    /// Return the (normalized) axis direction vector
    const VectorDouble& getAxisDirection() const { return _axisDir; }

    /// Dump a space in a string
    String toStringIdx(const AStringFormat* strfmt, Id idx = -1) const override;

    /// Return true if the given space is equal to me
    bool isEqual(const ASpace* space) const override;

  protected:
    /// Move the given space point by the given vector
    void _move(SpacePoint& p1, const VectorDouble& vec) const override;

    /// Return the geodesic distance between two space points (unrolled cylinder)
    double
      _getDistance(const SpacePoint& p1, const SpacePoint& p2, Id ispace = -1)
        const override;

    /// Return the distance between two space points with the given tensor
    double _getDistance(
      const SpacePoint& p1,
      const SpacePoint& p2,
      const Tensor& tensor,
      Id ispace = -1) const override;

    /// Return the distance in frequential domain between two space points with the given tensor
    double _getFrequentialDistance(
      const SpacePoint& p1,
      const SpacePoint& p2,
      const Tensor& tensor,
      Id ispace = -1) const override;

    /// Return the increment vector between two space points
    VectorDouble
      _getIncrement(const SpacePoint& p1, const SpacePoint& p2, Id ispace = -1)
        const override;

    /// Return the increment vector between two space points in a given vector
    void _getIncrementInPlace(
      const SpacePoint& p1,
      const SpacePoint& p2,
      VectorDouble& ptemp,
      Id ispace = -1) const override;

  private:
    /// Decompose a point into (height along axis, radial unit vector)
    void _decompose(const SpacePoint& p, double& h, VectorDouble& runit) const;

  private:
    /// Point through which the cylinder axis passes
    VectorDouble _axisOrigin;
    /// Normalized axis direction vector
    VectorDouble _axisDir;
    /// Cylinder radius
    double _radius;
  };
} // namespace gstlrn
