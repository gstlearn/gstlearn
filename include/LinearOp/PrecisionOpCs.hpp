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
#include "LinearOp/Cholesky.hpp"

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
                Model* model,
                int icov = 0,
                bool flagDecompose = false,
                const CGParam params = CGParam(),
                bool verbose = false);
  virtual ~PrecisionOpCs();

  // Interface for PrecisionOp class
  void evalDirect(const VectorDouble &vecin, VectorDouble &vecout) override;
  void evalSimulate(VectorDouble& whitenoise, VectorDouble& vecout) override;
  void evalInverse(VectorDouble& vecin, VectorDouble& vecout) override;
  void makeReady() override;

  double getLogDeterminant(int nbsimu = 1, int seed = 0) override;

  void evalDeriv(const VectorDouble& inv, VectorDouble& outv,int iapex,int igparam,const EPowerPT& power) override;
  void evalDerivOptim(VectorDouble& outv,int iapex,int igparam, const EPowerPT& power) override;
  //void evalDerivPoly(const VectorDouble& inv, VectorDouble& outv,int iapex,int igparam) override;
  void gradYQX(const VectorDouble & X, const VectorDouble &Y,VectorDouble& result, const EPowerPT& power) override;
  void gradYQXOptim(const VectorDouble & X, const VectorDouble &Y,VectorDouble& result, const EPowerPT& power) override;

  MatrixSparse* getQ() const { return _Q; }
  NF_Triplet getQToTriplet(bool flag_from_1 = false) const;

private:
  void _buildQ(bool flagDecompose = false);

private:
  MatrixSparse* _Q;
};
