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

#include "Basic/VectorNumT.hpp"
#include "LinearOp/ALinearOp.hpp"
#include "LinearOp/APreconditioner.hpp"
#include "LinearOp/CholeskySparse.hpp"
#include "LinearOp/IProj.hpp"
#include "Polynomials/Chebychev.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <memory>
#include <vector>

namespace gstlrn
{

  class Chebychev;

  class GSTLEARN_EXPORT MultiGridSolver: public APreconditioner
  {

  public:
    MultiGridSolver(
      double ratiochebmin = 0.1,
      double ratiochebmax = 1.1,
      size_t smoothIter = 4);

    MultiGridSolver(const MultiGridSolver& other) = delete;
    MultiGridSolver(MultiGridSolver&& other) noexcept = default;
    MultiGridSolver& operator=(const MultiGridSolver& other) = delete;
    MultiGridSolver& operator=(MultiGridSolver&& other) noexcept = default;
    virtual ~MultiGridSolver() = default;

    MultiGridSolver& compute(const ALinearOp& /*default*/) { return *this; }

    void prepare(size_t niterPower = 20);

    VectorDouble getMaxEigenValues() const { return _maxEigenValues; }

    const ALinearOp& getOperator(size_t lvl) const
    {
      return *_operators[lvl].get();
    }

    const std::vector<std::unique_ptr<const ALinearOp>>& getOperators() const
    {
      return _operators;
    }

    size_t getNLevels() const { return _nLevels; }

    void display() const;

#ifndef SWIG
    template<typename Rhs>
    Rhs solve(const Rhs& b) const
    {
      // Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
      // _vCycle(0, b, x);
      return b;
    }

    MultiGridSolver& analyzePattern(const ALinearOp& /*default*/)
    {
      return *this;
    }

    MultiGridSolver& factorize(const ALinearOp& /*default*/) { return *this; }

    void addLevel(
      std::unique_ptr<const ALinearOp>&& op,
      std::unique_ptr<const IProj>&& transferOp);
    void setCoarseSolver(const MatrixSparse& mat);
#endif

    void applyLevel(constvect rhs, vect x, size_t lvl) const;
    VectorDouble
      vCycle(const VectorDouble& rhs, const VectorDouble& u, size_t lvl) const;

  private:
    void _vCycle(size_t lvl, constvect f, vect u) const;

    size_t _nLevels;
    double _ratiochebmin;
    double _ratiochebmax;
    size_t _smoothIter;
    std::vector<std::unique_ptr<const IProj>> _transferOps;
    std::vector<std::unique_ptr<const ALinearOp>> _operators;
    std::vector<Chebychev> _chebys;
    VectorDouble _lMax;

    std::unique_ptr<const CholeskySparse> _coarseSolver;
    VectorDouble _maxEigenValues;
    mutable VectorVectorDouble _work;
    mutable VectorVectorDouble _res;
    mutable VectorVectorDouble _err_storage;
    mutable VectorVectorDouble _rhs_storage;
  };

} // namespace gstlrn
