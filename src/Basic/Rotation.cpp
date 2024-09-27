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
#include "Geometry/GeometryHelper.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"

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
    if (! Rotation::isMatrixRotation(rotmat, true)) return 1;
    _rotMat = rotmat;
    GH::rotationGetAnglesInPlace(local, _angles);
    _directToInverse();
    _checkRotForIdentity();
  }
  return 0;
}

int Rotation::setMatrixDirectVec(const VectorDouble& rotmat)
{
  if (! rotmat.empty())
  {
    if ((int) rotmat.size() != _rotMat.size())
      my_throw ("The argument 'rotmat' does not have same dimension as 'this'");
    MatrixSquareGeneral local(_nDim);
    local.setValues(rotmat);
    if (! Rotation::isMatrixRotation(local, true)) return 1;
    _rotMat = local;
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
  this->rotateDirect(inv.getVector(), outv.getVector());
}

void Rotation::rotateDirect(const std::vector<double>& inv, std::vector<double>& outv) const
{
  if (!_flagRot)
    outv = inv;
  else
   _rotMat.prodMatVecInPlace(inv, outv, false);
}

void Rotation::rotateInverse(const VectorDouble& inv, VectorDouble& outv) const
{
  this->rotateInverse(inv.getVector(), outv.getVector());
}

void Rotation::rotateInverse(const std::vector<double>& inv, std::vector<double>& outv) const
{
  if (!_flagRot)
    outv = inv;
  else
    _rotInv.prodMatVecInPlace(inv, outv, false);
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
  _flagRot = (! _rotMat.isIdentity());
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
bool Rotation::isMatrixRotation(const MatrixSquareGeneral& rotmat, bool verbose)
{

  /* Check product of matrix by its transpose and compare to unity matrix */

  int neq = rotmat.getNRows();
  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
    {
      double prod = 0.;
      for (int k = 0; k < neq; k++)
        prod += rotmat.getValue(i,k) * rotmat.getValue(j,k);
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
