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

#include "Basic/VectorNumT.hpp"
#include "Basic/WarningMacro.hpp"
#include "LinearOp/ACholesky.hpp"

#ifndef SWIG
#  include <Eigen/Core>
#  include <Eigen/Dense>
#  include <Eigen/src/Core/Matrix.h>
#endif

#ifndef SWIG
DISABLE_WARNING_PUSH
DISABLE_WARNING_COND_EXPR_CONSTANT
DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE
DISABLE_WARNING_DECLARATION_HIDE_GLOBAL
#  include <Eigen/Sparse>
DISABLE_WARNING_POP
#endif

namespace gstlrn
{
class MatrixSparse;
using Sp = Eigen::SparseMatrix<double>;

class GSTLEARN_EXPORT CholeskySparse: public ACholesky
{
public:
  CholeskySparse(const MatrixSparse& mat);
  CholeskySparse(const CholeskySparse& m);
  CholeskySparse& operator=(const CholeskySparse& m);
  virtual ~CholeskySparse();

  Id setMatrix(const MatrixSparse& mat);
  Id stdev(VectorDouble& vcur,
            const MatrixSparse* proj,
            bool flagStDev = false) const;

  double computeLogDeterminant() const override;
  Id addSolveX(const constvect vecin, vect vecout) const override;
  Id addInvLtX(const constvect vecin, vect vecout) const override;
  Id addLtX(const constvect vecin, vect vecout) const override;
  Id addLX(const constvect vecin, vect vecout) const override;
  Id addInvLX(const constvect vecin, vect vecout) const override;

private:
  void _clean();
  Id _prepare(const MatrixSparse& mat) const;
  Id _stdev(VectorDouble& vcur, const MatrixSparse* proj) const;

private:
  mutable Eigen::SimplicialLDLT<Sp>* _factor;
};
} // namespace gstlrn