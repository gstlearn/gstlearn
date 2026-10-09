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
#include "Space/SpaceCyl.hpp"
#include "Basic/AException.hpp"
#include "Basic/Tensor.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpacePoint.hpp"

#include <algorithm>
#include <cmath>
#include <memory>

namespace gstlrn
{
  namespace
  {
    double dot3(const VectorDouble& a, const VectorDouble& b)
    {
      return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    double norm3(const VectorDouble& a)
    {
      return sqrt(dot3(a, a));
    }
  } // namespace

  SpaceCyl::SpaceCyl(
    size_t ndim,
    const VectorDouble& axisOrigin,
    const VectorDouble& axisDir,
    double radius)
    : ASpace(ndim)
    , _axisOrigin(axisOrigin)
    , _axisDir(axisDir)
    , _radius(radius)
  {
    if (ndim != 3)
      my_throw("CYL is only implemented for ndim=3 (cylinder embedded in R3)");
    if (_axisOrigin.size() != 3) my_throw("axisOrigin must have 3 coordinates");
    if (_axisDir.size() != 3) my_throw("axisDir must have 3 coordinates");
    if (_radius <= 0.) my_throw("radius must be strictly positive");

    double n = norm3(_axisDir);
    if (n <= 0.) my_throw("axisDir must not be the null vector");
    for (Id i = 0; i < 3; i++) _axisDir[i] /= n;
  }

  SpaceCyl::SpaceCyl(const SpaceCyl& r)
    : ASpace(r)
    , _axisOrigin(r._axisOrigin)
    , _axisDir(r._axisDir)
    , _radius(r._radius)
  {
  }

  SpaceCyl& SpaceCyl::operator=(const SpaceCyl& r)
  {
    if (this != &r)
    {
      ASpace::operator=(r);
      _axisOrigin = r._axisOrigin;
      _axisDir = r._axisDir;
      _radius = r._radius;
    }
    return *this;
  }

  SpaceCyl::~SpaceCyl() {}

  ASpaceSharedPtr SpaceCyl::create(
    Id ndim,
    const VectorDouble& axisOrigin,
    const VectorDouble& axisDir,
    double radius)
  {
    return std::shared_ptr<SpaceCyl>(
      new SpaceCyl(ndim, axisOrigin, axisDir, radius));
  }

  String SpaceCyl::toStringIdx(const AStringFormat* strfmt, Id idx) const
  {
    std::stringstream sstr;
    sstr << ASpace::toStringIdx(strfmt, idx);
    if (strfmt == nullptr || strfmt->getLevel() == 1)
    {
      String suffix = (idx < 0) ? "" : (" [" + std::to_string(idx) + "]");
      sstr << "Cylinder Radius " << suffix << " = " << _radius << std::endl;
      sstr << "Axis Origin     " << suffix << " = (" << _axisOrigin[0] << ", "
           << _axisOrigin[1] << ", " << _axisOrigin[2] << ")" << std::endl;
      sstr << "Axis Direction  " << suffix << " = (" << _axisDir[0] << ", "
           << _axisDir[1] << ", " << _axisDir[2] << ")" << std::endl;
    }
    return sstr.str();
  }

  bool SpaceCyl::isEqual(const ASpace* space) const
  {
    if (!ASpace::isEqual(space)) return false;
    const auto* s = dynamic_cast<const SpaceCyl*>(space);
    return s != nullptr && _radius == s->_radius
        && _axisOrigin == s->_axisOrigin && _axisDir == s->_axisDir;
  }

  /// Decompose a point p into its height h along the axis and its
  /// (unit) radial direction runit, orthogonal to the axis.
  void SpaceCyl::_decompose(const SpacePoint& p, double& h, VectorDouble& runit)
    const
  {
    auto offset = static_cast<Id>(getOffset());

    VectorDouble v(3);
    for (Id i = 0; i < 3; i++) v[i] = p.getCoord(offset + i) - _axisOrigin[i];

    h = dot3(v, _axisDir);

    runit.resize(3);
    for (Id i = 0; i < 3; i++) runit[i] = v[i] - h * _axisDir[i];

    double n = norm3(runit);
    if (n > 0.)
      for (Id i = 0; i < 3; i++) runit[i] /= n;
  }

  void SpaceCyl::_move(SpacePoint& p1, const VectorDouble& vec) const
  {
    /// Note: moving a point by a raw cartesian vector does not, in
    /// general, keep it on the cylinder surface. It is the caller's
    /// responsibility to provide a displacement that is tangent to
    /// the surface if the constraint must be preserved.
    auto offset = static_cast<Id>(getOffset());
    auto ndim = static_cast<Id>(getNDim());
    for (Id i = offset; i < ndim + offset; i++)
    {
      p1.setCoord(i, p1.getCoord(i) + vec[i]);
    }
  }

  /**
   * Return the geodesic distance between two points lying on the
   * cylinder surface. Since a cylinder is a developable surface
   * (zero Gaussian curvature), unrolling it flattens the geodesics
   * into straight lines:
   *   d = sqrt( (radius * dtheta)^2 + dh^2 )
   * where dtheta is the angle between the two radial directions and
   * dh is the difference of heights along the axis.
   */
  double SpaceCyl::_getDistance(
    const SpacePoint& p1,
    const SpacePoint& p2,
    Id ispace) const
  {
    DECLARE_UNUSED(ispace)

    double h1, h2;
    VectorDouble r1, r2;
    _decompose(p1, h1, r1);
    _decompose(p2, h2, r2);

    double costheta = dot3(r1, r2);
    costheta = std::max(-1., std::min(1., costheta));
    double dtheta = acos(costheta);

    double dh = h2 - h1;
    double darc = _radius * dtheta;
    return sqrt(darc * darc + dh * dh);
  }

  double SpaceCyl::_getDistance(
    const SpacePoint& p1,
    const SpacePoint& p2,
    const Tensor& tensor,
    Id ispace) const
  {
    DECLARE_UNUSED(ispace)
    /// TODO : SpaceCyl::_getDistance with tensor (anisotropy on the
    /// unrolled (arc, height) plane is not yet implemented)
    DECLARE_UNUSED(tensor);
    return _getDistance(p1, p2, ispace);
  }

  double SpaceCyl::_getFrequentialDistance(
    const SpacePoint& p1,
    const SpacePoint& p2,
    const Tensor& tensor,
    Id ispace) const
  {
    DECLARE_UNUSED(ispace)
    /// TODO : SpaceCyl::_getFrequentialDistance
    DECLARE_UNUSED(p1);
    DECLARE_UNUSED(p2);
    DECLARE_UNUSED(tensor);
    return 0.;
  }

  VectorDouble SpaceCyl::_getIncrement(
    const SpacePoint& p1,
    const SpacePoint& p2,
    Id ispace) const
  {
    DECLARE_UNUSED(ispace)
    _getIncrementInPlace(p1, p2, _work1);
    return _work1;
  }

  /// Increment expressed in the raw cartesian coordinates
  /// (same convention as SpaceSN::_getIncrementInPlace)
  void SpaceCyl::_getIncrementInPlace(
    const SpacePoint& p1,
    const SpacePoint& p2,
    VectorDouble& ptemp,
    Id ispace) const
  {
    DECLARE_UNUSED(ispace)
    /// TODO : SpaceCyl::_getIncrementInPlace (raw cartesian difference,
    /// not a true geodesic transport on the surface)
    Id j = 0;
    auto offset = static_cast<Id>(getOffset());
    auto ndim = static_cast<Id>(getNDim());
    for (Id i = offset; i < ndim + offset; i++)
      ptemp[j++] = p2.getCoord(i) - p1.getCoord(i);
  }
} // namespace gstlrn
