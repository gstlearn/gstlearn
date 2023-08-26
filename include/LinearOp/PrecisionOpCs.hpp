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
#include "Matrix/LinkMatrixSparse.hpp"

class AMesh;
class ShiftOpCs;
class CovAniso;
class Model;

/** This class is just a specialization of PrecisionOp when the shift
* Operator is built with sparse (cs) matrices and therefore algebra can be performed with Cholesky.
* It allows to return the precision matrix as a cs. */

class GSTLEARN_EXPORT PrecisionOpCs : public PrecisionOp
{
public:
  PrecisionOpCs(ShiftOpCs* shiftop = nullptr,
                const CovAniso* cova = nullptr,
                bool verbose = false);
  PrecisionOpCs(const AMesh* mesh,
                Model* model,
                int icov = 0,
                bool verbose = false);
  virtual ~PrecisionOpCs();

  // Interface for PrecisionOp class
  void eval(const VectorDouble &inv, VectorDouble &outv) override;
  void simulateOneInPlace(VectorDouble& whitenoise, VectorDouble& result) override;
  void evalInvVect(VectorDouble& in, VectorDouble& result) override;
  double computeLogDet(int nbsimu = 1, int seed = 0) override;

  void evalDeriv(const VectorDouble& inv, VectorDouble& outv,int iapex,int igparam,const EPowerPT& power) override;
  void evalDerivOptim(VectorDouble& outv,int iapex,int igparam, const EPowerPT& power) override;
  //void evalDerivPoly(const VectorDouble& inv, VectorDouble& outv,int iapex,int igparam) override;
  void gradYQX(const VectorDouble & X, const VectorDouble &Y,VectorDouble& result, const EPowerPT& power) override;
  void gradYQXOptim(const VectorDouble & X, const VectorDouble &Y,VectorDouble& result, const EPowerPT& power) override;

#ifndef SWIG
  const cs* getQ() const { return _Q; }
#endif
  Triplet getQToTriplet(bool flag_from_1 = false) const;

private:
  void _buildQ();

private:
  cs* _Q;
  Cholesky _qChol;
};
