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
#include "Basic/Tensor.hpp"

#include "geoslib_f.h"
#include "geoslib_f_private.h"

#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"

Tensor::Tensor(unsigned int ndim)
:  AStringable(),
   _nDim(ndim),
   _tensorDirect(),
   _tensorInverse(),
   _radius(),
   _rotation(),
   _isotrop(true)
{
  init(ndim);
}

Tensor::Tensor(const Tensor& r)
:  AStringable(r),
   _nDim(r._nDim),
   _tensorDirect(r._tensorDirect),
   _tensorInverse(r._tensorInverse),
   _radius(r._radius),
   _rotation(r._rotation),
   _isotrop(r._isotrop)
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
    _radius = r._radius;
    _rotation = r._rotation;
    _isotrop = r._isotrop;
  }
  return *this;
}

Tensor::~Tensor()
{
}

void Tensor::init(int ndim)
{
  _radius.resize(_nDim, 1.);
  _rotation.resetFromSpaceDimension(_nDim);
  _rotation.setIdentity();
  _tensorDirect  = _rotation.getMatrixDirect();
  _tensorInverse = _rotation.getMatrixDirect();
  _isotrop = true;
}

String Tensor::toString(int /*level*/) const
{
  std::stringstream sstr;
  sstr << "Radius     = " << ut_vector_string(_radius) << std::endl;
  if (!_rotation.isIdentity())
    sstr << _rotation.toString() << std::endl;
  return sstr.str();
}

void Tensor::setRadius(double radius)
{
  if (ABS(radius) < EPSILON10)
    my_throw ("Ellipsoid radius cannot be null");
  ut_vector_fill(_radius, radius, static_cast<int> (_radius.size()));
  _isotrop = true;
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
  _updateIsotrop();
  _fillTensors();
}

void Tensor::setRadiusDir(unsigned int idim, double radius)
{
  if (idim >= _nDim)
    my_throw ("Wrong index of dimension");
  if (ABS(radius) < EPSILON10)
    my_throw ("Radius cannot be null");
  _radius[idim] = radius;
  _updateIsotrop();
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

VectorDouble Tensor::applyDirect(const VectorDouble& vec) const
{
  VectorDouble out = vec;
  _tensorDirect.prodVector(vec, out);
  return out;
}

VectorDouble Tensor::applyInverse(const VectorDouble& vec) const
{
  VectorDouble out = vec;
  _tensorInverse.prodVector(vec, out);
  return out;
}

void Tensor::_updateIsotrop()
{
  double rad0 = _radius[0];
  // for (const auto& rad1 : _radius) // (rad1 is forbidden under windows)
  for (const auto& r : _radius)
  {
    if (ABS(r - rad0) > EPSILON10 * (ABS(r) + ABS(rad0)))
    {
      _isotrop = false;
      return;
    }
  }
  _isotrop = true;
}

void Tensor::_fillTensors()
{
  _tensorDirect = _rotation.getMatrixInverse();
  _tensorDirect.multiplyRow(_radius);

  _tensorInverse = _rotation.getMatrixDirect();
  _tensorInverse.divideColumn(_radius);
}
