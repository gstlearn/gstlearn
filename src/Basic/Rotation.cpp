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
#include "Geometry/Rotation.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Matrix/MatrixSquare.hpp"

namespace gstlrn
{
Rotation::Rotation(size_t ndim)
  : AStringable()
  , _nDim(ndim)
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

Rotation& Rotation::operator=(const Rotation& r)
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

void Rotation::resetFromSpaceDimension(size_t ndim)
{
  _nDim    = ndim;
  _flagRot = false;
  _angles.resize(_nDim, 0.);
  _rotMat.reset(static_cast<Id>(_nDim), static_cast<Id>(_nDim));
  _rotMat.setIdentity();
  _rotInv.reset(static_cast<Id>(_nDim), static_cast<Id>(_nDim));
  _rotInv.setIdentity();
}

Id Rotation::setMatrixDirect(const MatrixSquare& rotmat)
{
  if (!rotmat.empty())
  {
    if (!_rotMat.isSameSize(rotmat))
      my_throw("The argument 'rotmat' does not have same dimension as 'this'");
    VectorDouble local = rotmat.getValues();
    if (!Rotation::isMatrixRotation(rotmat, true)) return 1;
    _rotMat = rotmat;
    GH::rotationGetAnglesInPlace(local, _angles);
    _directToInverse();
    _checkRotForIdentity();
  }
  return 0;
}

Id Rotation::setMatrixDirectVec(const VectorDouble& rotmat)
{
  if (!rotmat.empty())
  {
    if (static_cast<Id>(rotmat.size()) != _rotMat.size())
      my_throw("The argument 'rotmat' does not have same dimension as 'this'");
    MatrixSquare local(static_cast<Id>(_nDim));
    local.setValues(rotmat);
    if (!Rotation::isMatrixRotation(local, true)) return 1;
    _rotMat = local;
    GH::rotationGetAnglesInPlace(static_cast<Id>(_nDim), rotmat.data(), _angles.data());
    _directToInverse();
    _checkRotForIdentity();
  }
  return 0;
}

Id Rotation::setAngles(const VectorDouble& angles)
{
  if (!angles.empty())
  {
    if (angles.size() > _nDim)
      my_throw("Wrong dimension number for 'angles' argument");

    _angles = angles;
    _angles.resize(_nDim, 0.);
    if (_nDim == 2) _angles[1] = 0.;

    _local.resize(_nDim * _nDim);
    GH::rotationMatrixInPlace(static_cast<Id>(_nDim), _angles, _local);
    _rotMat.setValues(_local);
    _directToInverse();
    _checkRotForIdentity();
  }
  return 0;
}

Id Rotation::getDerivativesInPlace(std::vector<MatrixSquare>& res) const
{

  GH::rotationMatrixDerivativesInPlace(static_cast<Id>(_nDim), _angles, res);
  for (auto& dR: res)
  {
    dR.prodScalar(GV_PI / 180);
  }
  return 0;
}

std::vector<MatrixSquare> Rotation::getDerivatives() const
{
  std::vector<MatrixSquare> res;
  if (_nDim == 2)
  {
    res.resize(1);
    res[0].reset(2, 2);
  }
  else if (_nDim == 3)
  {
    res.resize(3);
    for (Id i = 0; i < 3; i++)
      res[i].reset(3, 3);
  }
  getDerivativesInPlace(res);
  return res;
}

void Rotation::setIdentity()
{
  for (Id idim = 0; idim < static_cast<Id>(_nDim); idim++)
    VH::fill(_angles, 0.);
  _rotMat.setIdentity();
  _rotInv.setIdentity();
  _checkRotForIdentity();
}

String Rotation::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (!_flagRot) return sstr.str();
  sstr << toVector("Rotation Angles        = ", _angles);

  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 0)
  {
    sstr << toMatrix("Direct Rotation Matrix", VectorString(), VectorString(),
                     true, static_cast<Id>(_nDim), static_cast<Id>(_nDim), _rotMat.getValues());
    sstr << toMatrix("Inverse Rotation Matrix", VectorString(), VectorString(),
                     true, static_cast<Id>(_nDim), static_cast<Id>(_nDim), _rotInv.getValues());
  }
  return sstr.str();
}

void Rotation::rotateDirect(const VectorDouble& inv, VectorDouble& outv) const
{
  if (!_flagRot)
    outv = inv;
  else
  {
    // Using the constvect interface allows resizing the input and output vectors
    // on the fly (avoiding copies)
    constvect cinv(inv.data(), _nDim);
    vect coutv(outv.data(), _nDim);
    _rotMat.prodMatVecInPlaceC(cinv, coutv, false);
  }
}

void Rotation::rotateInverse(const VectorDouble& inv, VectorDouble& outv) const
{
  if (!_flagRot)
    outv = inv;
  else
  {
    // Using the constvect interface allows resizing the input and output vectors
    // on the fly (avoiding copies)
    constvect cinv(inv.data(), _nDim);
    vect coutv(outv.data(), _nDim);
    _rotInv.prodMatVecInPlaceC(cinv, coutv, false);
  }
}

void Rotation::_recopy(const Rotation& r)
{
  _nDim    = r._nDim;
  _flagRot = r._flagRot;
  _angles  = r._angles;
  _rotMat  = r._rotMat;
  _rotInv  = r._rotInv;
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
  _flagRot = (!_rotMat.isIdentity());
}

bool Rotation::isSame(const Rotation& rot) const
{
  /* Find the minimum space dimension */

  Id ndim = static_cast<Id>(MIN(_nDim, rot.getNDim()));

  /* Compare the rotations */

  if (_flagRot != isRotated()) return 0;
  if (_flagRot)
  {
    for (Id idim = 0; idim < ndim; idim++)
      if (_angles[idim] != getAngle(idim)) return 0;
  }
  return 1;
}

/****************************************************************************/
/*!
 **  Check if a matrix is a rotation matrix
 **
 ** \return  true if the matrix is a rotation matrix; false otherwise
 **
 ** \param[in]  rotmat   Square matrix to be checked
 ** \param[in]  verbose  1 for the verbose option
 **
 ** \remark  A rotation matrix must be orthogonal with determinant equal to 1
 **
 *****************************************************************************/
bool Rotation::isMatrixRotation(const MatrixSquare& rotmat, bool verbose)
{

  /* Check product of matrix by its transpose and compare to unity matrix */

  auto neq = rotmat.getNRows();
  for (Id i = 0; i < neq; i++)
    for (Id j = 0; j < neq; j++)
    {
      double prod = 0.;
      for (Id k = 0; k < neq; k++)
        prod += rotmat.getValue(i, k) * rotmat.getValue(j, k);
      double comp = (i == j) ? 1 : 0.;
      if (ABS(prod - comp) > EPSILON6)
      {
        if (verbose)
          messerr("The element (A*At)[%d,%d] = %lf (should be %lf)", i + 1,
                  j + 1, prod, comp);
        return false;
      }
    }

  double deter = rotmat.determinant();
  if (ABS(deter - 1.) > EPSILON6)
  {
    if (verbose) messerr("The Determinant = %f (should be 1)", deter);
    return false;
  }
  return true;
}
} // namespace gstlrn