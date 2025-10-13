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

#include "API/SPDEParam.hpp"
#include "LinearOp/ASimulable.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "gstlearn_export.hpp"

#include <memory>
namespace gstlrn
{
class CholeskyDense;
class Db;
class Model;
class MatrixSymmetric;
/**
 * invNuggetOp
 *
 * This class is used to represent the inverse covariance matrix operator of a nugget effect.
 * It can be multivariate with anisotropy.
 * It stores the inverse submatrices for each location, allowing efficient computation of
 * the logdeterminant.
 * It inherits from ASimulable, allowing it to be used in simulations.
 */
class GSTLEARN_EXPORT InvNuggetOp: public virtual ASimulable, public virtual MatrixSparse
{
public:
  InvNuggetOp(const Db* dbin          = nullptr,
              Model* model            = nullptr,
              const SPDEParam& params = SPDEParam(),
              bool flagEigVals        = false);
  InvNuggetOp(const InvNuggetOp& m)            = delete;
  InvNuggetOp& operator=(const InvNuggetOp& m) = delete;
  virtual ~InvNuggetOp();
  Id getSize() const override;
  const MatrixSparse* getInvNuggetMatrix() const { return dynamic_cast<const MatrixSparse*>(this); }
  const MatrixSparse* cloneInvNuggetMatrix() const;
  const MatrixSparse* cloneCholNuggetMatrix() const;
  double computeLogDet(Id nMC = 1) const override;
  double getMinEigenValue() const { return _rangeEigenVal.first; }
  double getMaxEigenValue() const { return _rangeEigenVal.second; }
  std::pair<double, double> getRangeEigenValue() const { return _rangeEigenVal; }

protected:
  Id _addSimulateToDest(const constvect whitenoise, vect outv) const override;

private:
  void _updateMatrix(MatrixSymmetric& invsill,
                     CholeskyDense& cholsill,
                     Id ndef,
                     const VectorInt& position,
                     const VectorInt& identity);
  double _updateQuantities(MatrixSymmetric& sillsinv);
  void _buildInvNugget(const Db* dbin = nullptr, Model* model = nullptr, const SPDEParam& params = SPDEParam());

private:
  std::shared_ptr<MatrixSparse> _cholNuggetMatrix; // The Cholesky decomposition of the nugget matrix
  double _logDeterminant;                          // The log determinant of the inverse nugget matrix
  std::pair<double, double> _rangeEigenVal;        // The range of eigenvalues for the inverse nugget matrix
  bool _flagEigVals;                               // Flag to indicate if eigenvalues should be computed
};

#ifndef SWIG
// Function to build the inverse nugget matrix
GSTLEARN_EXPORT std::shared_ptr<MatrixSparse> buildInvNugget(Db* dbin, Model* model, const SPDEParam& params = SPDEParam());
#endif
} // namespace gstlrn
