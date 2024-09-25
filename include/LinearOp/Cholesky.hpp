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
  #include <Eigen/Core>
  #include <Eigen/Dense>
  #include <Eigen/src/Core/Matrix.h>
#endif

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

/// TODO : Inherit from ALinearOpEigenCG and remove evalInverse ?
class GSTLEARN_EXPORT Cholesky: public ALinearOp
{
public:
  Cholesky(const MatrixSparse* mat);
  Cholesky(const Cholesky &m) = delete;
  Cholesky& operator=(const Cholesky &m) = delete;
  virtual ~Cholesky();

  int getSize() const override;
  void evalInverse(const VectorDouble& vecin, VectorDouble& vecout) const;

  bool isValid() const { return _matCS != nullptr; }

  int  solve(const VectorDouble& b, VectorDouble& x) const;
  int  simulate(const VectorDouble& b, VectorDouble& x) const;
  #ifndef SWIG
  int solve(const constvect b, std::vector<double>& x) const;
  int solve(const constvect b, vect x) const;
  int simulate(const constvect b, vect x) const;
  int addSimulateToDest(const constvect b, vect x) const;
#endif
  int  stdev(VectorDouble& vcur, bool flagStDev = false) const;
  double getLogDeterminant() const;

#ifndef SWIG
protected:
  int _addToDest(const constvect inv, vect outv) const override;

private:
  void _clean();
  void _compute();

private:
  css *_S; // Cholesky decomposition (for Old-style Csparse storage)
  csn *_N; // Cholesky decomposition (for Old-style Csparse storage)
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > _cholSolver; // (for Eigen library storage)
  
  const MatrixSparse* _matCS; // Stored by compliance with ALinearOp. Not to be deleted
#endif
};
