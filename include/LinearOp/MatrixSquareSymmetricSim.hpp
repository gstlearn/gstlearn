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

#include "Basic/VectorNumT.hpp"
#include "LinearOp/ASimulable.hpp"
#include "Matrix/AMatrix.hpp"
#include <Eigen/Sparse>
class AMatrix;
class Cholesky;

/**
 * Square Symmetric matrices
 */
class GSTLEARN_EXPORT MatrixSquareSymmetricSim : public ASimulable
{
public:
  MatrixSquareSymmetricSim();
  MatrixSquareSymmetricSim(const MatrixSquareSymmetricSim &m) = delete;
  MatrixSquareSymmetricSim& operator=(const MatrixSquareSymmetricSim &m) = delete;
  MatrixSquareSymmetricSim(const AMatrix* m,bool inverse = true);
  virtual ~MatrixSquareSymmetricSim();

  const AMatrix* getMatrix() const {return _mat;}
  int getSize() const override { return _mat->getNRows();}
  bool isSparse() const {return _sparse;}
  bool isInverse() const {return _inverse;}
  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  bool isEmpty() const { return _empty;}


protected:
  virtual int _addSimulateToDest(const Eigen::VectorXd& whitenoise,
                                       Eigen::VectorXd& outv) const override;
  virtual int _addToDest(const Eigen::VectorXd& inv,
                               Eigen::VectorXd& outv) const override;


private:
  void _prepare() const; // modify only mutable objects and can be called from const method
  void _clear();
  
  int  _simulateSparse(const Eigen::VectorXd& whitenoise,
                             Eigen::VectorXd& outv) const;
private:
  const   AMatrix* _mat;
  bool             _inverse;
  bool             _empty;
  bool             _sparse;
 
  //quantities for Cholesky of Dense Matrix : TOTO : create ACholesky  for Cholesky
  // (renamed into CholeskySparse) and a new CholeskyDense class with a factory 
  //Choleskycompute(AMatrix*). ACholesky need an empty state (bool) 
  mutable Eigen::LLT<Eigen::MatrixXd>* _factorDense; // Cholesky decomposition (Eigen format)
  mutable VectorDouble _tl; // Lower triangular matrix (after Cholesky decomposition)
  mutable VectorDouble _xl; // Lower triangular matrix (inverse of _tl)
  ////////////////////////////////////////////////////////////////////////////
  mutable Cholesky* _factorSparse; // Cholesky decomposition
  
};