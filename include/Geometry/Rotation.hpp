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
#pragma once

#include "gstlearn_export.hpp"

#include "geoslib_define.h"
#include "Basic/AStringable.hpp"

#include "Matrix/MatrixSquareGeneral.hpp"

class GSTLEARN_EXPORT Rotation: public AStringable /// TODO : public ASpaceObject
{
public:
  Rotation(unsigned int ndim = 2);
  Rotation(const Rotation& r);
  Rotation& operator=(const Rotation& r);
  virtual ~Rotation();

  static bool isMatrixRotation(const MatrixSquareGeneral& rotmat, bool verbose);

  unsigned int getNDim() const { return _nDim; }
  bool isRotated() const { return _flagRot; }
  const MatrixSquareGeneral& getMatrixDirect() const { return _rotMat; }
  const MatrixSquareGeneral& getMatrixInverse() const { return _rotInv; }
  const VectorDouble& getAngles() const { return _angles; }
  double getAngle(int idim) const { return _angles[idim]; }

  void resetFromSpaceDimension(unsigned int ndim);
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;
  int setMatrixDirect(const MatrixSquareGeneral& rotmat);
  int setMatrixDirectVec(const VectorDouble& rotmat);
  int setAngles(const VectorDouble& angles);
  void setIdentity();
  void rotateDirect(const VectorDouble& inv, VectorDouble& outv) const;
  void rotateInverse(const VectorDouble& inv, VectorDouble& outv) const;
#ifndef SWIG
  void rotateDirect(const std::vector<double>& inv, std::vector<double>& outv) const;
  void rotateInverse(const std::vector<double>& inv, std::vector<double>& outv) const;
#endif // SWIG
  bool isIdentity() const { return !_flagRot; }
  bool isSame(const Rotation& rot) const;

  VectorDouble getMatrixDirectVec() const { return _rotMat.getValues(); }
  VectorDouble getMatrixInverseVec() const { return _rotInv.getValues(); }

  double getMatrixDirect(int idim, int jdim)  const { return _rotMat.getValue(idim, jdim); }
  double getMatrixInverse(int idim, int jdim) const { return _rotInv.getValue(idim, jdim); }

private:
  void _recopy(const Rotation& r);
  void _checkRotForIdentity();
  void _directToInverse();
  void _inverseToDirect();

private:
  unsigned int _nDim;
  bool _flagRot; // true if a Rotation is defined other than Identity
  VectorDouble    _angles;
  MatrixSquareGeneral _rotMat;
  MatrixSquareGeneral _rotInv;
};
