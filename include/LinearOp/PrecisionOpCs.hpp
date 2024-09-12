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

#include "LinearOp/Cholesky.hpp"
#include "gstlearn_export.hpp"
#include "LinearOp/PrecisionOp.hpp"

#ifndef SWIG
  #include <Eigen/src/Core/Matrix.h>
#endif

class AMesh;
class Cholesky;
class ShiftOpCs;
class CovAniso;
class Model;

/** This class is just a specialization of PrecisionOp when the shift
* Operator is built with sparse matrices and therefore algebra can be performed with Cholesky.
* It allows to return the precision matrix as a Sparse Matrix. */

class GSTLEARN_EXPORT PrecisionOpCs : public PrecisionOp
{
public:
  PrecisionOpCs(ShiftOpCs* shiftop = nullptr,
                const CovAniso* cova = nullptr,
                bool verbose = false);
  PrecisionOpCs(const AMesh* mesh,
                CovAniso* cova,
                bool verbose = false);
  virtual ~PrecisionOpCs();

  // Interface for PrecisionOp class
  #ifndef SWIG
  void evalInverse(const Eigen::VectorXd& vecin, Eigen::VectorXd& vecout) override;
  int _addSimulateToDest(const Eigen::VectorXd &whitenoise, Eigen::VectorXd& outv) const override;
  int _addToDest(const Eigen::VectorXd &inv, Eigen::VectorXd& outv) const override;
  #endif

  double getLogDeterminant(int nbsimu = 1, int seed = 0) override;
  
  //void evalDerivPoly(const VectorDouble& inv, VectorDouble& outv,int iapex,int igparam) override;
  #ifndef SWIG
  void evalDeriv(const Eigen::VectorXd& inv, Eigen::VectorXd& outv,int iapex,int igparam,const EPowerPT& power) override;
  void evalDerivOptim(Eigen::VectorXd& outv,int iapex,int igparam, const EPowerPT& power) override;
  void gradYQX(const Eigen::VectorXd & X, 
               const Eigen::VectorXd &Y,
               Eigen::VectorXd& result, const EPowerPT& power) override;
  void gradYQXOptim(const Eigen::VectorXd & X, const Eigen::VectorXd &Y,Eigen::VectorXd& result, const EPowerPT& power) override;
  #endif
  const MatrixSparse* getQ() const { return _Q; }

private:
  void _buildQ();
  MatrixSparse* _build_Q();

private:
  MatrixSparse* _Q;
  mutable Cholesky* _chol;
};
