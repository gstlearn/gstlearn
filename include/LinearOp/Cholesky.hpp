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

#include "Basic/WarningMacro.hpp"
#include "LinearOp/ALinearOp.hpp"
#include "Basic/VectorNumT.hpp"

#ifndef SWIG
DISABLE_WARNING_PUSH
DISABLE_WARNING_COND_EXPR_CONSTANT
DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE
DISABLE_WARNING_DECLARATION_HIDE_GLOBAL
#include <Eigen/Sparse>
DISABLE_WARNING_POP
#endif

class css; /// TODO : Dependency to csparse to be removed
class csn;
class MatrixSparse;

class GSTLEARN_EXPORT Cholesky: public ALinearOp
{
public:
  Cholesky(const MatrixSparse* mat);
  Cholesky(const Cholesky &m) = delete;
  Cholesky& operator=(const Cholesky &m) = delete;
  virtual ~Cholesky();

  int getSize() const override;
  void evalInverse(const VectorDouble& vecin, VectorDouble& vecout) const override;

  bool isValid() const { return _matCS != nullptr; }

  int  solve(const VectorDouble& b, VectorDouble& x) const;
  int  simulate(const VectorDouble& b, VectorDouble& x) const;
  int  stdev(VectorDouble& vcur, bool flagStDev = false) const;
  double getLogDeterminant() const;

protected:
  void _evalDirect(const VectorDouble& inv, VectorDouble& outv) const override;

private:
  void _clean();
  void _compute();

private:
#ifndef SWIG
  css *_S; // Cholesky decomposition (for Old-style Csparse storage)
  csn *_N; // Cholesky decomposition (for Old-style Csparse storage)
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > _cholSolver; // (for Eigen library storage)
#endif
  const MatrixSparse* _matCS; // Stored by compliance with ALinearOp. Not to be deleted
};
