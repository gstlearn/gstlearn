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
#include "geoslib_f_private.h"

#include "Matrix/MatrixSquareDiagonal.hpp"
#include "Basic/Tensor.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"

Tensor::Tensor(unsigned int ndim)
:  AStringable(),
   _nDim(ndim),
   _tensorDirect(),
   _tensorDirect2(),
   _tensorInverse(),
   _radius(),
   _rotation(),
   _isotropic(true)
{
  init(ndim);
}

Tensor::Tensor(const Tensor& r)
:  AStringable(r),
   _nDim(r._nDim),
   _tensorDirect(r._tensorDirect),
   _tensorDirect2(r._tensorDirect2),
   _tensorInverse(r._tensorInverse),
   _radius(r._radius),
   _rotation(r._rotation),
   _isotropic(r._isotropic)
{
}

Tensor& Tensor::operator=(const Tensor &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _nDim = r._nDim;
    _tensorDirect = r._tensorDirect;
    _tensorDirect2 = r._tensorDirect2;
    _tensorInverse = r._tensorInverse;
    _radius = r._radius;
    _rotation = r._rotation;
    _isotropic = r._isotropic;
  }
  return *this;
}

Tensor::~Tensor()
{
}

void Tensor::init(int ndim)
{
  _nDim = ndim;
  _radius.resize(_nDim, 1.);
  _rotation.resetFromSpaceDimension(_nDim);
  _rotation.setIdentity();
  _tensorDirect  = _rotation.getMatrixDirect();
  _tensorDirect2 = _rotation.getMatrixDirect();
  _tensorInverse = _rotation.getMatrixDirect();
  _isotropic = true;
}

String Tensor::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << "Radius     = " << VH::toString(_radius) << std::endl;
  if (!_rotation.isIdentity())
    sstr << _rotation.toString() << std::endl;
  return sstr.str();
}

void Tensor::setRadius(double radius)
{
  if (ABS(radius) < EPSILON10)
    my_throw ("Ellipsoid radius cannot be null");
  VH::fill(_radius, radius, static_cast<int> (_radius.size()));
  _isotropic = true;
  _fillTensors();
}

void Tensor::setRadiusVec(const VectorDouble& radius)
{
  if (radius.size() != _nDim || radius.size() == 0)
    my_throw ("Wrong dimension number for argument 'radius'");
  for (const auto& r : radius)
  {
    if (ABS(r) < EPSILON20)
      my_throw ("Radius cannot be null");
  }
  _radius = radius;
  _updateIsotropic();
  _fillTensors();
}

void Tensor::setRadiusDir(unsigned int idim, double radius)
{
  if (idim >= _nDim)
    my_throw ("Wrong index of dimension");
  if (ABS(radius) < EPSILON10)
    my_throw ("Radius cannot be null");
  _radius[idim] = radius;
  _updateIsotropic();
  _fillTensors();
}

void Tensor::setRotation(const Rotation& rot)
{
  if (rot.getNDim() != _nDim)
    my_throw ("Wrong dimension number of argument 'Rotation'");
  _rotation = rot;
  _fillTensors();
}

void Tensor::setRotationAngles(const VectorDouble& angles)
{
  /// TODO : Rotation angles in 2D
  if (_nDim > 2 && angles.size() != _nDim)
    my_throw("Dimension of argument 'angles' should match the Space Dimension");
  _rotation.setAngles(angles);
  _fillTensors();
}

void Tensor::setRotationAngle(unsigned int idim, double angle)
{
  if ((_nDim == 2 && idim != 0) ||
      (_nDim >  2 && idim >= _nDim))
    my_throw ("Wrong rank for Angle");
  VectorDouble angles = _rotation.getAngles();
  angles[idim] = angle;
  _rotation.setAngles(angles);
  _fillTensors();
}

VectorDouble Tensor::applyDirect(const VectorDouble& vec, int mode) const
{
  VectorDouble out = vec;
  if (mode == 1)
    _tensorDirect.prodVector(vec, out);
  else
    _tensorDirect2.prodVector(vec, out);
  return out;
}

void Tensor::applyInverseInPlace(const double* vec,double* out) const
{
    _tensorInverse.prodVector(vec, out);
}

void Tensor::applyInverseInPlace(const VectorDouble& vec,VectorDouble& out) const
{
    _tensorInverse.prodVector(vec, out);
}


VectorDouble Tensor::applyInverse(const VectorDouble& vec, int mode) const
{
  VectorDouble out = vec;
  if (mode == 1)
    _tensorInverse.prodVector(vec, out);
  else
    _tensorInverse.prodVector(vec, out);
  return out;
}

void Tensor::_updateIsotropic()
{
  double rad0 = _radius[0];
  // for (const auto& rad1 : _radius) // (rad1 is forbidden under windows)
  for (const auto& r : _radius)
  {
    if (ABS(r - rad0) > EPSILON10 * (ABS(r) + ABS(rad0)))
    {
      _isotropic = false;
      return;
    }
  }
  _isotropic = true;
}

void Tensor::_fillTensors()
{
  // Tensor = Radius %*% Rotation
  _tensorDirect = _rotation.getMatrixInverse();
  _tensorDirect.multiplyRow(_radius);

  // Tensor = Radius %*% Rotation (correct product)
  _tensorDirect2 = _tensorDirect;
  _tensorDirect2.transposeInPlace();

  // Tensor = Rotation %*% 1/Radius
  _tensorInverse = _rotation.getMatrixDirect();
  _tensorInverse.divideColumn(_radius);

  // Tensor = Rotation %*% 1/Radius (correct product)
  _tensorInverse2 = _tensorInverse;
  _tensorInverse2.transposeInPlace();
}
