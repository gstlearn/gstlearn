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
#include "LinearOp/ACholesky.hpp"
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

class GSTLEARN_EXPORT CholeskySparse: public ACholesky
{
public:
  CholeskySparse(const MatrixSparse* mat);
  CholeskySparse(const CholeskySparse& m)            = delete;
  CholeskySparse& operator=(const CholeskySparse& m) = delete;
  virtual ~CholeskySparse();

  int stdev(VectorDouble& vcur, bool flagStDev = false) const;

  double computeLogDeterminant() const override;
  int addSolveX(const constvect vecin, vect vecout) const override;
  int addInvLtX(const constvect vecin, vect vecout) const override;
  int addLtX(const constvect vecin, vect vecout) const override;
  int addLX(const constvect vecin, vect vecout) const override;
  int addInvLX(const constvect vecin, vect vecout) const override;

private:
  void _clean();
  void _prepare() const;

private:
  bool _flagEigen;

  // Old-style storage
  mutable css *_S; // Cholesky decomposition (for Old-style Csparse storage)
  mutable csn* _N; // Cholesky decomposition (for Old-style Csparse storage)

  // Eigen storage
  mutable Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > *_factor; // (Eigen library)
};
