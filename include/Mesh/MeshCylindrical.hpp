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
#include "Mesh/MeshEStandard.hpp"

namespace gstlrn
{
  /**
   * @brief Meshing of a cylindrical surface, built by triangulating the
   *        cylinder "unrolled" into (theta, z) coordinates.
   *
   * The cylinder is defined by:
   *  - a point through which the axis passes (_axisOrigin), default (0,0,0)
   *  - a unit direction vector for the axis (_axisDir), default (1,0,0)
   *  - a radius (_radius)
   *
   * Apices are stored parametrically as (theta [degrees], z), like
   * MeshSpherical stores (longitude, latitude). The embedded 3-D cartesian
   * coordinates are computed on demand via getEmbeddedCoorPerApex /
   * getEmbeddedCoorPerMesh.
   */
  class GSTLEARN_EXPORT MeshCylindrical: public MeshEStandard
  {
  public:
    MeshCylindrical();
    MeshCylindrical(const MeshCylindrical& m);
    MeshCylindrical& operator=(const MeshCylindrical& m);
    virtual ~MeshCylindrical();

    /// Interface to AStringable
    String toString(const AStringFormat* strfmt = nullptr) const override;

    /// ASerializable interface
    String getNFName() const override { return "MeshCylindrical"; }

    /// Interface to AMesh
    Id getVariety() const override { return 1; }

    Id getEmbeddedNDim() const override { return 3; }

    void getEmbeddedCoorPerMesh(Id imesh, Id ic, VectorDouble& coords)
      const override;
    void getEmbeddedCoorPerApex(Id iapex, VectorDouble& coords) const override;
    void getBarycenterInPlace(Id imesh, vect coord) const override;
    double getMeshSize(Id imesh) const override;

    static MeshCylindrical* create(
      const VectorDouble& vecTheta,
      const VectorDouble& vecZ,
      const VectorDouble& axisOrigin = {0., 0., 0.},
      const VectorDouble& axisDir = {1., 0., 0.},
      double radius = 1.,
      bool verbose = false);

    /**
     * Build the meshing from the list of theta (degrees) and z values.
     * theta=0 and theta=360 are automatically added and stitched together.
     */
    Id resetFromCylinder(
      const VectorDouble& vecTheta,
      const VectorDouble& vecZ,
      const VectorDouble& axisOrigin = {0., 0., 0.},
      const VectorDouble& axisDir = {1., 0., 0.},
      double radius = 1.,
      bool verbose = false);

    double getRadius() const { return _radius; }

    const VectorDouble& getAxisOrigin() const { return _axisOrigin; }

    const VectorDouble& getAxisDirection() const { return _axisDir; }

  private:
    void _computeFrame();
    void _toCartesian(double thetaDeg, double z, VectorDouble& xyz) const;

  private:
    VectorDouble _axisOrigin; // Point through which the axis passes
    VectorDouble _axisDir; // Normalized axis direction vector
    VectorDouble _uAxis; // Unit vector orthogonal to axis (theta = 0)
    VectorDouble _vAxis; // Unit vector orthogonal to axis and uAxis
    double _radius; // Cylinder radius
  };
} // namespace gstlrn
