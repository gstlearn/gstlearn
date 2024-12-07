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

#include "LinearOp/CholeskySparse.hpp"
#include "LinearOp/PrecisionOp.hpp"

class AMesh;
class ShiftOpMatrix;
class CovAniso;
class Model;

/** This class is just a specialization of PrecisionOp when the shift
* Operator is built with sparse matrices and therefore algebra can be performed with Cholesky.
* It allows to return the precision matrix as a Sparse Matrix. */

class GSTLEARN_EXPORT PrecisionOpMatrix : public PrecisionOp
{
public:
  PrecisionOpMatrix(ShiftOpMatrix* shiftop = nullptr,
                const CovAniso* cova = nullptr,
                bool verbose = false);
  PrecisionOpMatrix(const AMesh* mesh,
                CovAniso* cova,
                bool verbose = false);
  virtual ~PrecisionOpMatrix();

  // Interface for PrecisionOp class
#ifndef SWIG
  void evalInverse(const constvect vecin, std::vector<double>& vecout) override;
  int _addSimulateToDest(const constvect whitenoise, vect outv) const override;
  int _addToDest(const constvect inv, vect outv) const override;
#endif

  double getLogDeterminant(int nbsimu = 1) override;
  VectorDouble extractDiag() const override;

  //void evalDerivPoly(const VectorDouble& inv, VectorDouble& outv,int iapex,int igparam) override;
#ifndef SWIG
  void evalDeriv(const constvect inv,
                 vect outv,
                 int iapex,
                 int igparam,
                 const EPowerPT& power) override;
  void evalDerivOptim(vect outv,
                      int iapex,
                      int igparam,
                      const EPowerPT& power) override;
  void gradYQX(const constvect X,
               const constvect Y,
               vect result,
               const EPowerPT& power) override;
  void gradYQXOptim(const constvect X,
                    const constvect Y,
                    vect result,
                    const EPowerPT& power) override;
#endif
  const MatrixSparse* getQ() const { return _Q; }

private:
  void _buildQ();
  MatrixSparse* _build_Q();

private:
  MatrixSparse* _Q;
  mutable CholeskySparse* _chol;
};
