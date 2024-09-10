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
#include "LinearOp/PrecisionOp.hpp"

#ifndef SWIG
  #include <Eigen/src/Core/Matrix.h>
#endif

class AMesh;
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
                bool flagDecompose = false,
                bool verbose = false);
  PrecisionOpCs(const AMesh* mesh,
                CovAniso* cova,
                bool flagDecompose = false,
                bool verbose = false);
  virtual ~PrecisionOpCs();

  // Interface for PrecisionOp class
  #ifndef SWIG
  void evalSimulate(const Eigen::VectorXd& whitenoise, Eigen::VectorXd& vecout) override;
  void evalInverse(const Eigen::VectorXd& vecin, Eigen::VectorXd& vecout) override;
  #endif

  //TODO : required to call the method of the mother class from python?????
  void evalSimulate(const VectorDouble& in, VectorDouble &out) 
  {
    PrecisionOp::evalSimulate(in,out);
  }
  void makeReady() override;

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
  void _buildQ(bool flagDecompose = false);

private:
  MatrixSparse* _Q;
};
