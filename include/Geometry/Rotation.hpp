/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "GeometryHelper.hpp"

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

  unsigned int getNDim() const { return _nDim; }
  bool isRotated() const { return _flagRot; }
  const MatrixSquareGeneral& getMatrixDirect() const { return _rotMat; }
  const MatrixSquareGeneral& getMatrixInverse() const { return _rotInv; }
  const VectorDouble& getAngles() const { return _angles; }
  double getAngle(int idim) const { return _angles[idim]; }

  void resetFromSpaceDimension(unsigned int ndim);
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;
  int setMatrixDirect(const MatrixSquareGeneral& rotmat);
  int setMatrixDirectByVector(const VectorDouble& rotmat);
  int setMatrixDirectOldStyle(const double* rotmat);
  int setAngles(const VectorDouble& angles);
  void setIdentity();
  void rotateDirect(const VectorDouble& inv, VectorDouble& outv) const;
  void rotateInverse(const VectorDouble& inv, VectorDouble& outv) const;
  bool isIdentity() const { return !_flagRot; }
  bool isSame(const Rotation& rot) const;

  VectorDouble setDirection(int ndim,
                            const VectorDouble& angles = VectorDouble(),
                            double radius = 1.) const;

  VectorDouble getMatrixDirectByVector() const;
  VectorDouble getMatrixInverseByVector() const;

private:
  void _recopy(const Rotation& r);
  void _checkRot();
  void _directToInverse();
  void _inverseToDirect();

private:
  unsigned int _nDim;
  bool _flagRot; // true if a Rotation is defined other than Identity
  VectorDouble    _angles;
  MatrixSquareGeneral _rotMat;
  MatrixSquareGeneral _rotInv;
};
