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
#include "Mesh/MeshCylindrical.hpp"
#include "Basic/AException.hpp"
#include "geoslib_define.h"

#include <algorithm>
#include <cmath>

namespace gstlrn
{
  namespace
  {
    constexpr double PI_LOCAL = 3.14159265358979323846;

    double dot3(const VectorDouble& a, const VectorDouble& b)
    {
      return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
    }

    void cross3(const VectorDouble& a, const VectorDouble& b, VectorDouble& c)
    {
      c.resize(3);
      c[0] = a[1] * b[2] - a[2] * b[1];
      c[1] = a[2] * b[0] - a[0] * b[2];
      c[2] = a[0] * b[1] - a[1] * b[0];
    }

    double norm3(const VectorDouble& a)
    {
      return sqrt(dot3(a, a));
    }
  } // namespace

  MeshCylindrical::MeshCylindrical()
    : MeshEStandard()
    , _axisOrigin{0., 0., 0.}
    , _axisDir{1., 0., 0.}
    , _uAxis()
    , _vAxis()
    , _radius(1.)
  {
  }

  MeshCylindrical::MeshCylindrical(const MeshCylindrical& m)
    : MeshEStandard(m)
    , _axisOrigin(m._axisOrigin)
    , _axisDir(m._axisDir)
    , _uAxis(m._uAxis)
    , _vAxis(m._vAxis)
    , _radius(m._radius)
  {
  }

  MeshCylindrical& MeshCylindrical::operator=(const MeshCylindrical& m)
  {
    if (this != &m)
    {
      MeshEStandard::operator=(m);
      _axisOrigin = m._axisOrigin;
      _axisDir = m._axisDir;
      _uAxis = m._uAxis;
      _vAxis = m._vAxis;
      _radius = m._radius;
    }
    return *this;
  }

  MeshCylindrical::~MeshCylindrical() {}

  String MeshCylindrical::toString(const AStringFormat* strfmt) const
  {
    std::stringstream sstr;
    sstr << "Meshing of a Cylinder (unrolled into Theta-Z coordinates)"
         << "\n";
    sstr << "Radius         = " << _radius << "\n";
    sstr << "Axis Origin    = (" << _axisOrigin[0] << ", " << _axisOrigin[1]
         << ", " << _axisOrigin[2] << ")" << "\n";
    sstr << "Axis Direction = (" << _axisDir[0] << ", " << _axisDir[1] << ", "
         << _axisDir[2] << ")" << "\n";
    sstr << MeshEStandard::toString(strfmt);
    return sstr.str();
  }

  /// Build an orthonormal frame (_uAxis, _vAxis, _axisDir) such that
  /// theta = 0 points towards _uAxis and theta increases towards _vAxis
  void MeshCylindrical::_computeFrame()
  {
    double n = norm3(_axisDir);
    if (n <= 0.) my_throw("axisDir must not be the null vector");
    for (Id i = 0; i < 3; i++) _axisDir[i] /= n;

    // Pick an auxiliary vector not colinear with _axisDir
    VectorDouble tmp{0., 0., 1.};
    if (fabs(dot3(tmp, _axisDir)) > 0.9) tmp = VectorDouble{0., 1., 0.};

    cross3(_axisDir, tmp, _vAxis);
    double nv = norm3(_vAxis);
    for (Id i = 0; i < 3; i++) _vAxis[i] /= nv;

    cross3(_vAxis, _axisDir, _uAxis); // already unit norm
  }

  void
    MeshCylindrical::_toCartesian(double thetaDeg, double z, VectorDouble& xyz)
      const
  {
    double theta = thetaDeg * PI_LOCAL / 180.;
    double c = cos(theta);
    double s = sin(theta);
    xyz.resize(3);
    for (Id i = 0; i < 3; i++)
      xyz[i] = _axisOrigin[i] + z * _axisDir[i]
             + _radius * (c * _uAxis[i] + s * _vAxis[i]);
  }

  void MeshCylindrical::getEmbeddedCoorPerApex(Id iapex, VectorDouble& coords)
    const
  {
    double theta = getApexCoor(iapex, 0);
    double z = getApexCoor(iapex, 1);
    _toCartesian(theta, z, coords);
  }

  void MeshCylindrical::getEmbeddedCoorPerMesh(
    Id imesh,
    Id ic,
    VectorDouble& coords) const
  {
    double theta = getCoor(imesh, ic, 0);
    double z = getCoor(imesh, ic, 1);
    _toCartesian(theta, z, coords);
  }

  void MeshCylindrical::getBarycenterInPlace(Id imesh, vect coord) const
  {
    // Barycenter computed in the embedded 3-D space (more meaningful on a
    // curved surface than an average of theta,z, which would be biased
    // close to the seam).
    auto ncorner = getNApexPerMesh();
    VectorDouble sum(3, 0.);
    VectorDouble xyz(3);
    for (Id ic = 0; ic < ncorner; ic++)
    {
      getEmbeddedCoorPerMesh(imesh, ic, xyz);
      for (Id i = 0; i < 3; i++) sum[i] += xyz[i];
    }
    for (Id i = 0; i < 3; i++) coord[i] = sum[i] / static_cast<double>(ncorner);
  }

  /// Real (developed) triangle area, computed on the unrolled cylinder in
  /// the (arclength, z) plane, with arclength = radius * theta(radians)
  double MeshCylindrical::getMeshSize(Id imesh) const
  {
    double t0 = getCoor(imesh, 0, 0) * PI_LOCAL / 180. * _radius;
    double z0 = getCoor(imesh, 0, 1);
    double t1 = getCoor(imesh, 1, 0) * PI_LOCAL / 180. * _radius;
    double z1 = getCoor(imesh, 1, 1);
    double t2 = getCoor(imesh, 2, 0) * PI_LOCAL / 180. * _radius;
    double z2 = getCoor(imesh, 2, 1);
    return 0.5 * fabs(((t1 - t0) * (z2 - z0)) - ((t2 - t0) * (z1 - z0)));
  }

  MeshCylindrical* MeshCylindrical::create(
    const VectorDouble& vecTheta,
    const VectorDouble& vecZ,
    const VectorDouble& axisOrigin,
    const VectorDouble& axisDir,
    double radius,
    bool verbose)
  {
    auto* mesh = new MeshCylindrical();
    (void)mesh->resetFromCylinder(
      vecTheta, vecZ, axisOrigin, axisDir, radius, verbose);
    return mesh;
  }

  /**
   * Build the meshing of the cylinder, unrolled into (theta, z) coordinates.
   *
   * Follows the same logic as the R prototype 'build_triangu':
   *  1. theta=0 and theta=360 are forced into the list of angles so the
   *     mesh covers the full circumference.
   *  2. The (theta, z) grid is a full Cartesian product ; since a Delaunay
   *     triangulation of a full rectangular grid is simply obtained by
   *     splitting each cell into 2 triangles (as done for regular grids in
   *     meshes_turbo_2D_grid_build / MeshETurbo), we do exactly that,
   *     instead of calling a generic (and here unnecessary) Delaunay solver.
   *  3. theta=0 and theta=360 represent the same meridian on the cylinder:
   *     every triangle index pointing to theta=360 is redirected to the
   *     corresponding node at theta=0 (same z).
   *  4. The (now unused) theta=360 nodes are dropped and indices remapped.
   */
  Id MeshCylindrical::resetFromCylinder(
    const VectorDouble& vecThetaIn,
    const VectorDouble& vecZIn,
    const VectorDouble& axisOrigin,
    const VectorDouble& axisDir,
    double radius,
    bool verbose)
  {
    DECLARE_UNUSED(verbose)
    if (axisOrigin.size() != 3 || axisDir.size() != 3)
      my_throw("axisOrigin and axisDir must have 3 coordinates");
    if (radius <= 0.) my_throw("radius must be strictly positive");

    _axisOrigin = axisOrigin;
    _axisDir = axisDir;
    _radius = radius;
    _computeFrame();

    // Force presence of theta=0 and theta=360 (edges of the unrolled
    // cylinder), sort, and remove duplicates
    VectorDouble vecTheta = vecThetaIn;
    vecTheta.push_back(0.);
    vecTheta.push_back(360.);
    std::sort(vecTheta.begin(), vecTheta.end());
    vecTheta.erase(
      std::unique(vecTheta.begin(), vecTheta.end()), vecTheta.end());

    VectorDouble vecZ = vecZIn;
    std::sort(vecZ.begin(), vecZ.end());
    vecZ.erase(std::unique(vecZ.begin(), vecZ.end()), vecZ.end());

    auto nt = static_cast<Id>(vecTheta.size());
    auto nz = static_cast<Id>(vecZ.size());
    if (nt < 2 || nz < 2)
      my_throw(
        "At least 2 distinct Theta values (besides 0/360) and 2 Z values "
        "are required");

    // Full grid of (theta, z) combinations: node index = it * nz + iz
    auto nodeIndex = [nz](Id it, Id iz) { return (it * nz) + iz; };
    Id nnodes = nt * nz;

    VectorDouble apicesRaw(nnodes * 2);
    for (Id it = 0; it < nt; it++)
      for (Id iz = 0; iz < nz; iz++)
      {
        Id ip = nodeIndex(it, iz);
        apicesRaw[(ip * 2) + 0] = vecTheta[it];
        apicesRaw[(ip * 2) + 1] = vecZ[iz];
      }

    // Triangulate each grid cell into 2 triangles (equivalent, on a full
    // rectangular grid, to a Delaunay triangulation)
    Id ncorner = 3;
    Id ncell = (nt - 1) * (nz - 1);
    VectorInt meshesRaw(ncell * 2 * ncorner);
    Id imesh = 0;
    for (Id it = 0; it < nt - 1; it++)
      for (Id iz = 0; iz < nz - 1; iz++)
      {
        Id i00 = nodeIndex(it, iz);
        Id i10 = nodeIndex(it + 1, iz);
        Id i01 = nodeIndex(it, iz + 1);
        Id i11 = nodeIndex(it + 1, iz + 1);

        // Triangle 1: (i00, i10, i11)
        meshesRaw[(imesh * ncorner) + 0] = i00;
        meshesRaw[(imesh * ncorner) + 1] = i10;
        meshesRaw[(imesh * ncorner) + 2] = i11;
        imesh++;

        // Triangle 2: (i00, i11, i01)
        meshesRaw[(imesh * ncorner) + 0] = i00;
        meshesRaw[(imesh * ncorner) + 1] = i11;
        meshesRaw[(imesh * ncorner) + 2] = i01;
        imesh++;
      }

    // Stitching theta=0 and theta=360 together: every node at theta=360
    // is redirected to the corresponding node at theta=0 (same z)
    VectorInt associate(nnodes);
    for (Id it = 0; it < nt; it++)
      for (Id iz = 0; iz < nz; iz++)
      {
        Id ip = nodeIndex(it, iz);
        if (vecTheta[it] == 360.)
          associate[ip] = nodeIndex(0, iz); // theta=0, same z
        else
          associate[ip] = ip;
      }

    for (Id& idx: meshesRaw) idx = associate[idx];

    // Drop the now-redundant theta=360 nodes and remap indices
    VectorInt oldToNew(nnodes, -1);
    Id nNewNodes = 0;
    VectorDouble apicesFinal;
    for (Id it = 0; it < nt; it++)
    {
      if (vecTheta[it] == 360.) continue; // drop
      for (Id iz = 0; iz < nz; iz++)
      {
        Id ip = nodeIndex(it, iz);
        oldToNew[ip] = nNewNodes;
        apicesFinal.push_back(vecTheta[it]);
        apicesFinal.push_back(vecZ[iz]);
        nNewNodes++;
      }
    }

    for (Id& idx: meshesRaw) idx = oldToNew[idx];

    // Store into the generic MeshEStandard storage (theta, z parametric
    // coordinates are stored as if they were "flat" 2D coordinates;
    // the true embedded 3-D position is reconstructed on demand via
    // getEmbeddedCoorPerApex / getEmbeddedCoorPerMesh)
    return resetFromVectors(
      2, ncorner, apicesFinal, meshesRaw, /* byCol */ false, verbose);
  }
} // namespace gstlrn
