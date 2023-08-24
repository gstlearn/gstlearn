/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/AStringable.hpp"
#include "Geometry/Rotation.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"

class GSTLEARN_EXPORT Tensor : public AStringable/// TODO : public ASpaceObject
{
public:
  Tensor(unsigned int ndim = 2);
  Tensor(const Tensor &r);
  Tensor& operator= (const Tensor &r);
  virtual ~Tensor();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(int ndim);
  void setTensorDirect (const MatrixSquareGeneral& tensor) { _tensorDirect  = tensor; }
  void setTensorInverse(const MatrixSquareGeneral& tensor) { _tensorInverse = tensor; }

  void setRadius(double radius); // Make it isotropic
  void setRadiusVec(const VectorDouble& radius);
  void setRadiusDir(unsigned int idim, double radius);

  void setRotation(const Rotation& rot);
  void setRotationAngles(const VectorDouble& angles);
  void setRotationAngle(unsigned int idim, double angle);

  const VectorDouble&    getAngles()        const { return  _rotation.getAngles(); }
  const MatrixSquareGeneral& getTensorDirect()  const { return  _tensorDirect; }
  const MatrixSquareGeneral& getTensorInverse() const { return  _tensorInverse; }
  const VectorDouble&    getRadius()        const { return  _radius; }
  const Rotation&        getRotation()      const { return  _rotation; }
  const MatrixSquareGeneral& getMatrixDirect()  const { return  _rotation.getMatrixDirect(); }
  const MatrixSquareGeneral& getMatrixInverse() const { return  _rotation.getMatrixInverse(); }
  bool                   isIsotropic()        const { return  _isotropic; }
  bool                   hasRotation()      const { return !_rotation.isIdentity(); }

  VectorDouble applyDirect (const VectorDouble& vec) const;
  VectorDouble applyInverse(const VectorDouble& vec) const;
  void applyInverseInPlace(const VectorDouble& vec, VectorDouble& out) const;
  void applyDirectInPlace(const VectorDouble &vec, VectorDouble &out) const;

private:
  void _updateIsotropic();
  void _fillTensors();

private:
  unsigned int        _nDim;     /// Number of dimensions
  MatrixSquareGeneral _tensorDirect; /// Direct Tensor matrix (definite positive)
  MatrixSquareGeneral _tensorInverse; /// Inverse Tensor matrix (definite positive)
  VectorDouble        _radius;   /// Ellipsoid radius
  Rotation            _rotation; /// Ellipsoid rotation
  bool                _isotropic;  /// True if the tensor is isotropic
};
