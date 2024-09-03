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

#include "Matrix/MatrixSparse.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"

#ifndef SWIG
  #include <Eigen/Core>
  #include <Eigen/Dense>
#endif

class Model;
/**
 * Class for the precision matrix of the latent field in SPDE (matricial form)
 */
class GSTLEARN_EXPORT PrecisionOpMultiMatrix :  public MatrixSparse ,  public PrecisionOpMulti
{
public:
  PrecisionOpMultiMatrix(Model* model = nullptr, 
                   const std::vector<const AMesh*>& meshes = std::vector<const AMesh*>());
  PrecisionOpMultiMatrix(const PrecisionOpMulti &m)= delete;
  PrecisionOpMultiMatrix& operator= (const PrecisionOpMulti &m)= delete;
  virtual ~PrecisionOpMultiMatrix();

  #ifndef SWIG
  
  protected:
  int    _addToDest(const Eigen::VectorXd& inv,
                        Eigen::VectorXd& outv) const override;
  #endif

  private:
  MatrixSparse _prepareMatrixNoStat(int icov, const MatrixSparse* Q);
  MatrixSparse _prepareMatrixStationary(int icov, const MatrixSparse* Q);
  void _prepareMatrix();
  void _makeReady() override;
  void _buildQop() override;

};
