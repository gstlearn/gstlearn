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

#include "Basic/Vector.hpp"
#include "Basic/Rotation.hpp"
#include "MatrixC/MatrixCSGeneral.hpp"
#include "Basic/AStringable.hpp"

class Tensor : public AStringable/// TODO : public ASpaceObject
{
public:
  Tensor(unsigned int ndim = 2);
  Tensor(const Tensor &r);
  Tensor& operator= (const Tensor &r);
  virtual ~Tensor();

  void init(int ndim);
  virtual String toString(int level = 0) const override;

  void setTensorDirect (const MatrixCSGeneral& tensor) { _tensorDirect  = tensor; }
  void setTensorInverse(const MatrixCSGeneral& tensor) { _tensorInverse = tensor; }

  void setRadius(double radius); // Make it isotropic
  void setRadius(const VectorDouble& radius);
  void setRadius(unsigned int idim, double radius);

  void setRotation(const Rotation& rot);
  void setRotationAngles(const VectorDouble& angles);
  void setRotationAngle(unsigned int idim, double angle);

  const VectorDouble&    getAngles()        const { return  _rotation.getAngles(); }
  const MatrixCSGeneral& getTensorDirect()  const { return  _tensorDirect; }
  const MatrixCSGeneral& getTensorInverse() const { return  _tensorInverse; }
  const VectorDouble&    getRadius()        const { return  _radius; }
  const Rotation&        getRotation()      const { return  _rotation; }
  const MatrixCSGeneral& getMatrixDirect()  const { return  _rotation.getMatrixDirect(); }
  const MatrixCSGeneral& getMatrixInverse() const { return  _rotation.getMatrixInverse(); }
  bool                   isIsotrop()        const { return  _isotrop; }
  bool                   hasRotation()      const { return !_rotation.isIdentity(); }

  VectorDouble applyDirect (const VectorDouble& vec) const;
  VectorDouble applyInverse(const VectorDouble& vec) const;

private:
  void _updateIsotrop();
  void _fillTensors();

private:
  unsigned int    _nDim;     /// Number of dimensions
  MatrixCSGeneral _tensorDirect; /// Direct Tensor matrix (definite positive)
  MatrixCSGeneral _tensorInverse; /// Inverse Tensor matrix (definite positive)
  VectorDouble    _radius;   /// Ellipsoid radius
  Rotation        _rotation; /// Ellipsoid rotation
  bool            _isotrop;  /// True if the tensor is isotropic
};

