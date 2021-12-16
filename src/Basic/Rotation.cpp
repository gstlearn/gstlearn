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
#include "Basic/Rotation.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

Rotation::Rotation(unsigned int ndim)
  : _nDim(ndim)
  , _flagRot(false)
  , _angles()
  , _rotMat()
  , _rotInv()
{
  resetFromSpaceDimension(ndim);
}

Rotation::Rotation(const Rotation& r)
{
  _recopy(r);
}

Rotation& Rotation::operator= (const Rotation& r)
{
  if (this != &r)
  {
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
  if (! rotmat.isEmpty())
  {
    if (! _rotMat.isSameSize(rotmat))
      my_throw ("The argument 'rotmat' does not have same dimension as 'this'");
    VectorDouble local = rotmat.getValues();
    if (! is_matrix_rotation(_nDim, local.data(), 1)) return 1;
    _rotMat.setValues(local);
    ut_angles_from_rotation_matrix(local.data(), _nDim, _angles.data());
    _directToInverse();
    _checkRot();
  }
  return 0;
}

int Rotation::setMatrixDirectByVector(const VectorDouble& rotmat)
{
  if (! rotmat.empty())
  {
    if ((int) rotmat.size() != _rotMat.getNTotal())
      my_throw ("The argument 'rotmat' does not have same dimension as 'this'");
    setMatrixDirectOldStyle(rotmat.data());
  }
  return 0;
}

VectorDouble Rotation::getMatrixDirectByVector() const
{
  return _rotMat.getValues();
}

VectorDouble Rotation::getMatrixInverseByVector() const
{
  return _rotMat.getValues();
}

int Rotation::setMatrixDirectOldStyle(const double* rotmat)
{
  if (rotmat != nullptr)
  {
    if (! is_matrix_rotation(_nDim, rotmat, 1)) return 1;
    _rotMat.setValues(rotmat);
    ut_angles_from_rotation_matrix(rotmat, _nDim, _angles.data());
    _directToInverse();
    _checkRot();
  }
  return 0;
}

int Rotation::setAngles(const VectorDouble& angles)
{
  if (! angles.empty())
  {
    if (angles.size() > _nDim)
      my_throw("Wrong dimension number for 'angles' argument");

    _angles.resize(_nDim,0.);
    int nval = MAX((int) _nDim, static_cast<int> (angles.size()));
    for (int idim = 0; idim < nval; idim++)
      _angles[idim] = angles[idim];
    if (_nDim == 2) _angles[1] = 0.;

    VectorDouble local = VectorDouble(_nDim * _nDim);
    ut_rotation_matrix(_nDim, _angles.data(), local.data());
    _rotMat.setValues(local);
    _directToInverse();
    _checkRot();
  }
  return 0;
}

void Rotation::setIdentity()
{
  for (int idim = 0; idim < (int) _nDim; idim++)
    ut_vector_fill(_angles,0.);
  _rotMat.setIdentity();
  _rotInv.setIdentity();
  _checkRot();
}

String Rotation::toString(int level) const
{
  std::stringstream sstr;
  if (!_flagRot) return sstr.str();
  sstr << toVector("Rotation Angles        = ",_angles);

  if (! IFFFF(level) && level > 0)
  {
    sstr << toMatrix("Direct Rotation Matrix", VectorString(), VectorString(),
                     false, _nDim, _nDim, _rotMat.getValues());
    sstr << toMatrix("Inverse Rotation Matrix", VectorString(), VectorString(),
                     false, _nDim, _nDim, _rotInv.getValues());
  }
  return sstr.str();
}

void Rotation::rotateDirect(const VectorDouble& in, VectorDouble& out) const
{
  if (!_flagRot)
    out = in;
  else
  _rotMat.prodVector(in, out);
}

void Rotation::rotateInverse(const VectorDouble& in, VectorDouble& out) const
{
  if (!_flagRot)
    out = in;
  else
    _rotInv.prodVector(in, out);
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
  if (_rotInv.invert())
  {
    messerr("Error in the inversion of the rotation matrix");
    messerr("The Rotation is cancelled");
  }
}

void Rotation::_inverseToDirect()
{
  _rotMat = _rotInv;
  if (_rotMat.invert())
  {
    messerr("Error in the inversion of the rotation matrix");
    messerr("The Rotation is cancelled");
  }
}

void Rotation::_checkRot()
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

VectorDouble Rotation::setDirection(int ndim,
                                    const VectorDouble& angles,
                                    double radius) const
{
  VectorDouble codir(ndim);
  (void) ut_angles_to_codir(ndim,1,angles,codir);
  if (radius != 1.)
    for (int idim = 0; idim < ndim; idim++) codir[idim] *= radius;
  return codir;
}
