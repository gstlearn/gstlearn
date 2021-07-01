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
#pragma once

#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"
#include "MatrixC/MatrixCSGeneral.hpp"

class Rotation: public AStringable /// TODO : public ASpaceObject
{
public:
  Rotation(unsigned int ndim = 2);
  Rotation(const Rotation& r);
  Rotation& operator=(const Rotation& r);
  virtual ~Rotation();

  unsigned int getNDim() const { return _nDim; }
  bool isRotated() const { return _flagRot; }
  const MatrixCSGeneral& getMatrixDirect() const { return _rotMat; }
  const MatrixCSGeneral& getMatrixInverse() const { return _rotInv; }
  const VectorDouble& getAngles() const { return _angles; }
  double getAngles(int idim) const { return _angles[idim]; }

  void init(unsigned int ndim);
  virtual String toString(int level = 0) const override;
  int setMatrixDirect(const MatrixCSGeneral& rotmat);
  int setMatrixDirect(const VectorDouble& rotmat);
  int setMatrixDirect(const double* rotmat);
  int setAngles(const VectorDouble& angles);
  void setIdentity();
  void rotateDirect(const VectorDouble& in, VectorDouble& out) const;
  void rotateInverse(const VectorDouble& in, VectorDouble& out) const;
  bool isIdentity() const { return !_flagRot; }
  bool isSame(const Rotation& rot) const;

  VectorDouble setDirection(int ndim,
                            const VectorDouble& angles = VectorDouble(),
                            double radius = 1.) const;

private:
  void _recopy(const Rotation& r);
  void _checkRot();
  void _directToInverse();
  void _inverseToDirect();

private:
  unsigned int _nDim;
  bool _flagRot; // true if a Rotation is defined other than Identity
  VectorDouble    _angles;
  MatrixCSGeneral _rotMat;
  MatrixCSGeneral _rotInv;
};
