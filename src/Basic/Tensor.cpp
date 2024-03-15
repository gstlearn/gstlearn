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
#include "geoslib_f_private.h"

#include "Basic/Tensor.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"

Tensor::Tensor(unsigned int ndim)
:  AStringable(),
   _nDim(ndim),
   _tensorDirect(),
   _tensorInverse(),
   _tensorDirect2(),
   _tensorInverse2(),
   _tensorDirectSwap(),
   _radius(),
   _rotation(),
   _isotropic(true),
   _flagDefinedBySquare(false)
{
  init(ndim);
}

Tensor::Tensor(const Tensor& r)
:  AStringable(r),
   _nDim(r._nDim),
   _tensorDirect(r._tensorDirect),
   _tensorInverse(r._tensorInverse),
   _tensorDirect2(r._tensorDirect2),
   _tensorInverse2(r._tensorInverse2),
   _tensorDirectSwap(r._tensorDirectSwap),
   _radius(r._radius),
   _rotation(r._rotation),
   _isotropic(r._isotropic),
   _flagDefinedBySquare(r._flagDefinedBySquare)
{
}

Tensor& Tensor::operator=(const Tensor &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _nDim = r._nDim;
    _tensorDirect = r._tensorDirect;
    _tensorInverse = r._tensorInverse;
    _tensorDirect2 = r._tensorDirect2;
    _tensorInverse2 = r._tensorInverse2;
    _tensorDirectSwap = r._tensorDirectSwap;
    _radius = r._radius;
    _rotation = r._rotation;
    _isotropic = r._isotropic;
    _flagDefinedBySquare = r._flagDefinedBySquare;
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
  _tensorDirect   = _rotation.getMatrixDirect();
  _tensorInverse  = _rotation.getMatrixInverse();
  // Squared tensor are easy to calculate as rotation is identity and radius is 1
  _tensorDirect2  = _rotation.getMatrixInverse();
  _tensorInverse2 = _rotation.getMatrixInverse();
  _tensorDirectSwap = _rotation.getMatrixDirect();
  _isotropic = true;
}

String Tensor::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << "Radius     = " << VH::toStringAsVD(_radius) << std::endl;
  if (!_rotation.isIdentity())
    sstr << _rotation.toString() << std::endl;
  return sstr.str();
}

void Tensor::setRadiusIsotropic(double radius)
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

/**
 * This functions defines jointly the rotation anisotropy and ranges.
 * It allows initiating the tensor only once (saves time)
 * @param angles Vector of rotation angles (optional)
 * @param radius Vector of ranges (optional)
 */
void Tensor::setRotationAnglesAndRadius(const VectorDouble& angles, const VectorDouble& radius)
{
  if (! angles.empty())
  {
    if (_nDim > 2 && angles.size() != _nDim)
      my_throw("Dimension of argument 'angles' should match the Space Dimension");
    _rotation.setAngles(angles);
  }

  if (!radius.empty())
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

  }
  _fillTensors();
}

VectorDouble Tensor::applyDirect(const VectorDouble& vec) const
{
  VectorDouble out = vec;
  _tensorDirect.prodMatVecInPlace(vec, out);
  return out;
}

void Tensor::applyInverseInPlace(const VectorDouble &vec, VectorDouble &out) const
{
  _tensorInverse.prodMatVecInPlace(vec, out);
}

void Tensor::applyInverse2InPlace(const VectorDouble &vec, VectorDouble &out) const
{
  _tensorInverse2.prodMatVecInPlace(vec, out);
}

void Tensor::applyDirectInPlace(const VectorDouble &vec, VectorDouble &out) const
{
  _tensorDirect.prodMatVecInPlace(vec, out);
}

void Tensor::applyDirectSwapInPlace(const VectorDouble &vec, VectorDouble &out) const
{
  _tensorDirectSwap.prodMatVecInPlace(vec, out);
}


VectorDouble Tensor::applyInverse(const VectorDouble& vec) const
{
  VectorDouble out = vec;
  _tensorInverse.prodMatVecInPlace(vec, out);
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
  _tensorDirect = _rotation.getMatrixDirect();
  _tensorDirect.multiplyColumn(_radius);

  // Tensor = t(Rotation) %*% 1/Radius
  _tensorInverse = _rotation.getMatrixInverse();
  _tensorInverse.divideRow(_radius);

  // Square of the Direct tensor
  _tensorDirect2 = MatrixSquareSymmetric(_nDim);
  _tensorDirect2.prodMatMatInPlace(&_tensorDirect, &_tensorDirect, false, true);

  _tensorDirectSwap = _rotation.getMatrixDirect();
  _tensorDirectSwap.multiplyRow(_radius);

  // Inverse of the Direct squared tensor
  _direct2ToInverse2();
}

void Tensor::_direct2ToInverse2()
{
  _tensorInverse2 = _tensorDirect2;
  _tensorInverse2.invert();
}

void Tensor::setTensorDirect2(const MatrixSquareSymmetric& tensor)
{
  _tensorDirect2 = tensor;
  _direct2ToInverse2();
  _flagDefinedBySquare = true;
}
