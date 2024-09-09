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

#include "Matrix/MatrixSparse.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"

#ifndef SWIG
  #include <Eigen/src/Core/Matrix.h>
  #include <Eigen/Core>
  #include <Eigen/Dense>
#endif

class Model;
/**
 * Class for the precision matrix of the latent field in SPDE (matricial form)
 */
class GSTLEARN_EXPORT PrecisionOpMultiMatrix : public PrecisionOpMulti
{
public:
  PrecisionOpMultiMatrix(Model* model               = nullptr,
                         const VectorMeshes& meshes = VectorMeshes());
  PrecisionOpMultiMatrix(const PrecisionOpMulti &m)= delete;
  PrecisionOpMultiMatrix& operator= (const PrecisionOpMulti &m)= delete;
  virtual ~PrecisionOpMultiMatrix();

  const MatrixSparse* getQ() const;
  private:
  #ifndef SWIG
  virtual int _addToDestImpl(const Eigen::VectorXd &vecin,Eigen::VectorXd &vecout) const override;
  #endif
  MatrixSparse _prepareMatrixNoStat(int icov, const MatrixSparse* Q) const;
  MatrixSparse _prepareMatrixStationary(int icov, const MatrixSparse* Q) const;
  void _prepareMatrix();
  void _makeReady() override;
  void _buildQop() override;
  bool _isSingle() const { return _getNVar() == 1 && _getNCov() == 1;}

private:
  MatrixSparse _Q;
};
