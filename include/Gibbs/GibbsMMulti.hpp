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

#include "Gibbs/GibbsMulti.hpp"
#include "LinearOp/CholeskySparse.hpp"

namespace gstlrn
{
class Db;
class Model;
class MatrixSparse;

class GSTLEARN_EXPORT GibbsMMulti: public GibbsMulti
{
public:
  GibbsMMulti();
  GibbsMMulti(Db* db, Model* model);
  GibbsMMulti(const GibbsMMulti& r);
  GibbsMMulti& operator=(const GibbsMMulti& r);
  virtual ~GibbsMMulti();

  void update(VectorVectorDouble& y, Id isimu, Id ipgs, Id iter) override;
  Id covmatAlloc(bool verbose, bool verboseTimer = false) override;

  void setEps(double eps) { _eps = eps; }
  void cleanup() override;

  void setFlagStoreInternal(bool flagStoreInternal);

private:
  Id _getNVar() const;
  Id _getSize() const;
  void _storeWeights(Id icol);
  void _storeWeightsMS(Id icol, NF_Triplet& NF_T);
  void _getWeights(Id icol) const;
  void _calculateWeights(Id icol,
                         VectorDouble& b,
                         double tol = EPSILON3) const;
  Id _storeAllWeights(bool verbose = false);
  Id _getSizeOfWeights(const VectorDouble& weights) const;
  double _getEstimate(Id ipgs, Id icol, const VectorVectorDouble& y) const;
  void _allocate();
  double _getVariance(Id icol) const;
  Id _getColumn(Id iact, Id ivar) const;
  void _splitCol(Id icol, Id* iact, Id* ivar) const;
  void _updateStatWeights(Id* nzero);

private:
  std::shared_ptr<MatrixSparse> _Cmat;
  CholeskySparse* _CmatChol;
  double _eps;
  bool _flagStoreInternal;

  VectorVectorDouble _areas;
  MatrixSparse* _matWgt;

  mutable VectorDouble _weights;
};
} // namespace gstlrn