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

  VectorDouble applyDirect (const VectorDouble& vec, int mode = 1) const;
  VectorDouble applyInverse(const VectorDouble& vec, int mode = 1) const;

private:
  void _updateIsotropic();
  void _fillTensors();

private:
  unsigned int    _nDim;     /// Number of dimensions
  MatrixSquareGeneral _tensorDirect; /// Direct Tensor matrix (definite positive)
  MatrixSquareGeneral _tensorDirect2; /// Direct Tensor matrix (definite positive)
  MatrixSquareGeneral _tensorInverse; /// Inverse Tensor matrix (definite positive)
  MatrixSquareGeneral _tensorInverse2; /// Inverse Tensor matrix (definite positive)
  VectorDouble    _radius;   /// Ellipsoid radius
  Rotation        _rotation; /// Ellipsoid rotation
  bool            _isotropic;  /// True if the tensor is isotropic
};
