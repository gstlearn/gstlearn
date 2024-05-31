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
#include "Geometry/BiTargetCheckDistance.hpp"
#include "geoslib_old_f.h"

#include "Basic/VectorHelper.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Space/SpaceTarget.hpp"

BiTargetCheckDistance::BiTargetCheckDistance(double radius,
                                             const VectorDouble coeffs,
                                             const VectorDouble angles)
  : ABiTargetCheck(),
    _ndim(2),
    _flagAniso(false),
    _flagRotation(false),
    _radius(TEST),
    _anisoCoeffs(),
    _anisoRotMat(),
    _dist(TEST),
    _movingIncr(),
    _movingAux()
{
  _radius = radius;
  _anisoCoeffs.resize(_ndim);
  _anisoRotMat.resize(_ndim * _ndim);
  if (! coeffs.empty())
  {
    _ndim = (int) coeffs.size();

    //    _flagAniso = (ut_vector_constant(coeffs)) ? 0 : 1;
    _flagAniso = true;
    _anisoCoeffs = coeffs;

    if (! angles.empty())
    {
      _flagRotation = (VH::isConstant(angles, 0.)) ? false : true;
      GH::rotationMatrixInPlace(_ndim, angles, _anisoRotMat);
    }
    else
    {
      GH::rotationMatrixIdentityInPlace(_ndim, _anisoRotMat);
    }
  }
  else
  {
    VH::fill(_anisoCoeffs, 1., _ndim);
    GH::rotationMatrixIdentityInPlace(_ndim, _anisoRotMat);
  }

  _movingIncr.resize(_ndim);
  _movingAux.resize(_ndim);
}

BiTargetCheckDistance::BiTargetCheckDistance(const BiTargetCheckDistance &r)
  : ABiTargetCheck(r),
    _ndim(r._ndim),
    _flagAniso(r._flagAniso),
    _flagRotation(r._flagRotation),
    _radius(r._radius),
    _anisoCoeffs(r._anisoCoeffs),
    _anisoRotMat(r._anisoRotMat),
    _dist(r._dist),
    _movingIncr(r._movingIncr),
    _movingAux(r._movingAux)
{
}

BiTargetCheckDistance& BiTargetCheckDistance::operator=(const BiTargetCheckDistance &r)
{
  if (this != &r)
  {
    ABiTargetCheck::operator=(r);
    _ndim = r._ndim;
    _flagAniso = r._flagAniso;
    _flagRotation = r._flagRotation;
    _radius = r._radius;
    _anisoCoeffs = r._anisoCoeffs;
    _anisoRotMat = r._anisoRotMat;
    _dist = r._dist;
    _movingIncr = r._movingIncr;
    _movingAux = r._movingAux;
  }
  return *this;
}

BiTargetCheckDistance::~BiTargetCheckDistance()
{
}

BiTargetCheckDistance* BiTargetCheckDistance::create(double radius,
                                                   const VectorDouble coeffs,
                                                   const VectorDouble angles)
{
  return new BiTargetCheckDistance(radius, coeffs, angles);
}

String BiTargetCheckDistance::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  if (!FFFF(_radius))
  {
    if (!_flagAniso)
    {
      sstr << "Maximum horizontal distance         = " << _radius
           << std::endl;
    }
    else
    {
      VectorDouble ranges(_ndim);
      for (int idim = 0; idim < _ndim; idim++)
        ranges[idim] = _radius * _anisoCoeffs[idim];
      sstr << toMatrix("Anisotropic Ranges :", VectorString(), VectorString(),
                      true, 1, _ndim, ranges);

      if (_flagRotation)
      {
        sstr << toMatrix("Anisotropy Rotation :", VectorString(),
                        VectorString(), true, _ndim, _ndim, _anisoRotMat);
      }
    }
  }
  return sstr.str();
}

double BiTargetCheckDistance::getNormalizedDistance(const VectorDouble& dd) const
{
  _movingIncr = dd;

  _calculateDistance();

  return _dist;
}

void BiTargetCheckDistance::_calculateDistance() const
{
  int ndim = getNDim();

  /* Anisotropic neighborhood */

  if (_flagAniso)
  {

    /* Rotated anisotropy ellipsoid */

    if (_flagRotation)
    {
      matrix_product_safe(1, ndim, ndim, _movingIncr.data(),
                          _anisoRotMat.data(), _movingAux.data());
      _movingIncr = _movingAux;
    }
    for (int idim = 0; idim < ndim; idim++)
      _movingIncr[idim] /= _anisoCoeffs[idim];
  }

  /* Calculate the distance */

  matrix_product_safe(1, ndim, 1, _movingIncr.data(), _movingIncr.data(), &_dist);
  _dist = sqrt(_dist);
}

bool BiTargetCheckDistance::isOK(const SpaceTarget &T1,
                                const SpaceTarget &T2) const
{
  int ndim = getNDim();
  for (int idim = 0; idim < ndim; idim++)
    _movingIncr[idim] = T1.getCoord(idim) - T2.getCoord(idim);

  _calculateDistance();

  return _dist <= _radius;
}
