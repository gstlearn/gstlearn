/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "Geometry/BiPointCheckDistance.hpp"
#include "Space/SpacePoint.hpp"

BiPointCheckDistance::BiPointCheckDistance(double radius,
                                           const VectorDouble coeffs,
                                           const VectorDouble angles)
    : ABiPointCheck(),
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
      GH::rotationInit(_ndim, angles.data(), _anisoRotMat.data());
    }
    else
    {
      GH::rotationIdentity(_ndim, _anisoRotMat.data());
    }
  }
  else
  {
    VH::fill(_anisoCoeffs, 1., _ndim);
    GH::rotationIdentity(_ndim, _anisoRotMat.data());
  }

  _movingIncr.resize(_ndim);
  _movingAux.resize(_ndim);
}

BiPointCheckDistance::BiPointCheckDistance(const BiPointCheckDistance &r)
    : ABiPointCheck(r),
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

BiPointCheckDistance& BiPointCheckDistance::operator=(const BiPointCheckDistance &r)
{
  if (this != &r)
  {
    ABiPointCheck::operator=(r);
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

BiPointCheckDistance::~BiPointCheckDistance()
{
}

BiPointCheckDistance* BiPointCheckDistance::create(double radius,
                                                   const VectorDouble coeffs,
                                                   const VectorDouble angles)
{
  return new BiPointCheckDistance(radius, coeffs, angles);
}

String BiPointCheckDistance::toString(const AStringFormat* /*strfmt*/) const
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
                      true, _ndim, 1, ranges);

      if (_flagRotation)
      {
        sstr << toMatrix("Anisotropy Rotation :", VectorString(),
                        VectorString(), true, _ndim, _ndim, _anisoRotMat);
      }
    }
  }
  return sstr.str();
}

double BiPointCheckDistance::getNormalizedDistance(const VectorDouble& dd) const
{
  _movingIncr = dd;

  _calculateDistance();

  return _dist;
}

void BiPointCheckDistance::_calculateDistance() const
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

bool BiPointCheckDistance::isOK(const SpacePoint &P1,
                                const SpacePoint &P2,
                                int iech1,
                                int iech2) const
{
  int ndim = getNDim();
  for (int idim = 0; idim < ndim; idim++)
    _movingIncr[idim] = P1.getCoord(idim) - P2.getCoord(idim);

  _calculateDistance();

  return _dist <= _radius;
}
