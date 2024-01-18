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
#include "geoslib_old_f.h"

#include "Geometry/Rotation.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"

Rotation::Rotation(unsigned int ndim)
  : AStringable(),
    _nDim(ndim)
  , _flagRot(false)
  , _angles()
  , _rotMat()
  , _rotInv()
{
  resetFromSpaceDimension(ndim);
}

Rotation::Rotation(const Rotation& r)
  : AStringable(r)
{
  _recopy(r);
}

Rotation& Rotation::operator= (const Rotation& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _recopy(r);
  }
  return *this;
}

Rotation::~Rotation()
{
}

void Rotation::resetFromSpaceDimension(unsigned int ndim)
{
  _nDim = ndim;
  _flagRot = false;
  _angles.resize(_nDim, 0.);
  _rotMat.reset(_nDim , _nDim);
  _rotMat.setIdentity();
  _rotInv.reset(_nDim , _nDim);
  _rotInv.setIdentity();
}

int Rotation::setMatrixDirect(const MatrixSquareGeneral& rotmat)
{
  if (! rotmat.empty())
  {
    if (! _rotMat.isSameSize(rotmat))
      my_throw ("The argument 'rotmat' does not have same dimension as 'this'");
    VectorDouble local = rotmat.getValues();
    if (! is_matrix_rotation(_nDim, local.data(), 1)) return 1;
    _rotMat = rotmat;
    GH::rotationGetAnglesInPlace(local, _angles);
    _directToInverse();
    _checkRotForIdentity();
  }
  return 0;
}

int Rotation::setMatrixDirectByVector(const VectorDouble& rotmat)
{
  if (! rotmat.empty())
  {
    if ((int) rotmat.size() != _rotMat.size())
      my_throw ("The argument 'rotmat' does not have same dimension as 'this'");
    if (! is_matrix_rotation(_nDim, rotmat.data(), 1)) return 1;
    _rotMat.setValues(rotmat);
    GH::rotationGetAnglesInPlace(_nDim, rotmat.data(), _angles.data());
    _directToInverse();
    _checkRotForIdentity();
  }
  return 0;
}

int Rotation::setAngles(const VectorDouble& angles)
{
  if (! angles.empty())
  {
    if (angles.size() > _nDim)
      my_throw("Wrong dimension number for 'angles' argument");

    _angles = angles;
    _angles.resize(_nDim,0.);
    if (_nDim == 2) _angles[1] = 0.;

    VectorDouble local = VectorDouble(_nDim * _nDim);
    GH::rotationMatrixInPlace(_nDim, _angles, local);
    _rotMat.setValues(local);
    _directToInverse();
    _checkRotForIdentity();
  }
  return 0;
}

void Rotation::setIdentity()
{
  for (int idim = 0; idim < (int) _nDim; idim++)
    VH::fill(_angles,0.);
  _rotMat.setIdentity();
  _rotInv.setIdentity();
  _checkRotForIdentity();
}

String Rotation::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (!_flagRot) return sstr.str();
  sstr << toVector("Rotation Angles        = ",_angles);

  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 0)
  {
    sstr << toMatrix("Direct Rotation Matrix", VectorString(), VectorString(),
                     true, _nDim, _nDim, _rotMat.getValues());
    sstr << toMatrix("Inverse Rotation Matrix", VectorString(), VectorString(),
                     true, _nDim, _nDim, _rotInv.getValues());
  }
  return sstr.str();
}

void Rotation::rotateDirect(const VectorDouble& inv, VectorDouble& outv) const
{
  if (!_flagRot)
    outv = inv;
  else
   _rotMat.prodVectorInPlace(inv, outv);
}

void Rotation::rotateInverse(const VectorDouble& inv, VectorDouble& outv) const
{
  if (!_flagRot)
    outv = inv;
  else
    _rotInv.prodVectorInPlace(inv, outv);
}

void Rotation::_recopy(const Rotation &r)
{
  _nDim = r._nDim;
  _flagRot = r._flagRot;
  _angles = r._angles;
  _rotMat = r._rotMat;
  _rotInv = r._rotInv;
}

void Rotation::_directToInverse()
{
  _rotInv = _rotMat;
  _rotInv.transposeInPlace();
}

void Rotation::_inverseToDirect()
{
  _rotMat = _rotInv;
  _rotMat.transposeInPlace();
}

void Rotation::_checkRotForIdentity()
{
  if (_rotMat.isIdentity())
    _flagRot = false;
  else
    _flagRot = true;
}

bool Rotation::isSame(const Rotation& rot) const
{
  /* Find the minimum space dimension */

  int ndim = MIN(_nDim, rot.getNDim());

  /* Compare the rotations */

  if (_flagRot != isRotated()) return 0;
  if (_flagRot)
  {
    for (int idim = 0; idim < ndim; idim++)
      if (_angles[idim] != getAngle(idim)) return 0;
  }
  return 1;
}
