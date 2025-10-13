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

#include "Basic/AStringable.hpp"
#include "geoslib_define.h"

#include "Matrix/MatrixSquare.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT Rotation: public AStringable /// TODO : public ASpaceObject
{
public:
  Rotation(size_t ndim = 2);
  Rotation(const Rotation& r);
  Rotation& operator=(const Rotation& r);
  virtual ~Rotation();

  static bool isMatrixRotation(const MatrixSquare& rotmat, bool verbose);

  size_t getNDim() const { return _nDim; }
  bool isRotated() const { return _flagRot; }
  const MatrixSquare& getMatrixDirect() const { return _rotMat; }
  const MatrixSquare& getMatrixInverse() const { return _rotInv; }
  const VectorDouble& getAngles() const { return _angles; }
  double getAngle(Id idim) const { return _angles[idim]; }
  Id getDerivativesInPlace(std::vector<MatrixSquare>& res) const;
  std::vector<MatrixSquare> getDerivatives() const;
  void resetFromSpaceDimension(size_t ndim);
  String toString(const AStringFormat* strfmt = nullptr) const override;
  Id setMatrixDirect(const MatrixSquare& rotmat);
  Id setMatrixDirectVec(const VectorDouble& rotmat);
  Id setAngles(const VectorDouble& angles);
  void setIdentity();
  void rotateDirect(const VectorDouble& inv, VectorDouble& outv) const;
  void rotateInverse(const VectorDouble& inv, VectorDouble& outv) const;
  bool isIdentity() const { return !_flagRot; }
  bool isSame(const Rotation& rot) const;

  VectorDouble getMatrixDirectVec() const { return _rotMat.getValues(); }
  VectorDouble getMatrixInverseVec() const { return _rotInv.getValues(); }

  double getMatrixDirect(Id idim, Id jdim) const { return _rotMat.getValue(idim, jdim); }
  double getMatrixInverse(Id idim, Id jdim) const { return _rotInv.getValue(idim, jdim); }

private:
  void _recopy(const Rotation& r);
  void _checkRotForIdentity();
  void _directToInverse();
  void _inverseToDirect();

private:
  size_t _nDim;
  bool _flagRot; // true if a Rotation is defined other than Identity
  VectorDouble _angles;
  MatrixSquare _rotMat;
  MatrixSquare _rotInv;
  mutable VectorDouble _local;
};
} // namespace gstlrn