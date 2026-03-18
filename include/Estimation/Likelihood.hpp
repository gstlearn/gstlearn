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

#include "Estimation/ALikelihood.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class Db;
class ModelGeneric;

class GSTLEARN_EXPORT Likelihood: public ALikelihood
{
public:
  Likelihood(ModelGeneric* model, const Db* db, bool reml = false);
  Likelihood(const Likelihood& r);
  Likelihood& operator=(const Likelihood& r);
  virtual ~Likelihood();

  static Likelihood* createForOptim(ModelGeneric* model,
                                    const Db* db,
                                    bool reml    = false,
                                    bool verbose = false);
  void evalGrad(vect res) override;
  MatrixSparse& getQMat() const override;


private:
  void _fillGradCovMat(RankHandler& rkh, const covmaptype& gradcov);
  void _updateModel(bool verbose = false) override;
  void _computeCm1X() override;
  void _computeCm1Yc() override;
  double _computeLogDet() const override;
  void _solveQ(constvect inv, vect outv) const override;
private:
  std::shared_ptr<MatrixSymmetric> _cov;
  CholeskyDense _covChol;
  MatrixSymmetric _gradCovMat;
  VectorDouble _temp;
  MatrixSquare _gradCovMatTimesInvCov;
};
} // namespace gstlrn
