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
#include "Geometry/Rotation.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixSymmetric.hpp"

class GSTLEARN_EXPORT Tensor : public AStringable /// TODO : public ASpaceObject
{
public:
  Tensor(unsigned int ndim = 2);
  Tensor(const Tensor &r);
  Tensor& operator= (const Tensor &r);
  virtual ~Tensor();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(int ndim);
  void setTensorDirect2(const MatrixSymmetric& tensor);

  void setRadiusIsotropic(double radius);
  void setRadiusVec(const VectorDouble& radius);
  void setRadiusDir(unsigned int idim, double radius);

  void setRotation(const Rotation& rot);
  void setRotationAngles(const VectorDouble& angles);
  void setRotationAngle(unsigned int idim, double angle);

  void setRotationAnglesAndRadius(const VectorDouble& angles = VectorDouble(),
                                  const VectorDouble& radius = VectorDouble());

  const VectorDouble&          getAngles()         const { return  _rotation.getAngles(); }
  const MatrixSquare&   getTensorDirect()   const { return  _tensorDirect; }
  const MatrixSquare&   getTensorInverse()  const { return  _tensorInverse; }
  const MatrixSymmetric& getTensorDirect2()  const { return  _tensorDirect2; }
  const VectorDouble&          getRadius()         const { return  _radius; }
  double getAngle(int i)         const { return  getAngles()[i]; }
  double getRadius(int i)         const { return  _radius[i]; }
  const Rotation&              getRotation()       const { return  _rotation; }
  const MatrixSquare&   getMatrixDirect()   const { return  _rotation.getMatrixDirect(); }
  const MatrixSquare&   getMatrixInverse()  const { return  _rotation.getMatrixInverse(); }
  bool                         isIsotropic()       const { return  _isotropic; }
  bool                         hasRotation()       const { return !_rotation.isIdentity(); }

  double getValue(int idim, int jdim) const { return _rotation.getMatrixDirect(idim, jdim); }

  VectorDouble applyDirect (const VectorDouble& vec) const;
  VectorDouble applyInverse(const VectorDouble& vec) const;
  void applyInverseInPlace(constvect vec, vect out) const;

  void applyInverseInPlace(const VectorDouble& vec, VectorDouble& out) const;
  void applyInverse2InPlace(const VectorDouble& vec, VectorDouble& out) const;
  void applyDirectInPlace(const VectorDouble &vec, VectorDouble &out) const;
  void applyDirectSwapInPlace(const VectorDouble &vec, VectorDouble &out) const;

  bool isFlagDefinedByInverse2() const { return _flagDefinedBySquare; }

private:
  void _updateIsotropic();
  void _fillTensors();
  void _direct2ToInverse2();

private:
  unsigned int _nDim;                      /// Number of dimensions
  MatrixSquare   _tensorDirect;     /// Direct Tensor matrix (definite positive)
  MatrixSquare   _tensorInverse;    /// Inverse Tensor matrix (definite positive)
  MatrixSymmetric _tensorDirect2;    /// Square of Direct tensor
  MatrixSymmetric _tensorInverse2;   /// Inverse of squared direct tensor
  MatrixSquare   _tensorDirectSwap; /// If tensor = Radius x Rot, tensorSwap = Rot x Radius
  VectorDouble _radius;                    /// Ellipsoid radius
  Rotation _rotation;                      /// Ellipsoid rotation
  bool _isotropic;                         /// True if the tensor is isotropic
  bool _flagDefinedBySquare;               /// True if Tensor has been defined using the squared elements (direct)
};
