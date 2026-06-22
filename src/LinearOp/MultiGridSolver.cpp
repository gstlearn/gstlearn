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

#include "LinearOp/MultiGridSolver.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "LinearOp/CholeskySparse.hpp"
#include "LinearOp/IProj.hpp"
#include "LinearOp/LinearOpHelper.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Polynomials/Chebychev.hpp"

namespace gstlrn
{

  MultiGridSolver::MultiGridSolver(
    double ratiochebmin,
    double ratiochebmax,
    size_t smoothIter)
    : _nLevels(1)
    , _ratiochebmin(ratiochebmin)
    , _ratiochebmax(ratiochebmax)
    , _smoothIter(smoothIter)
  {
  }

  void MultiGridSolver::addLevel(
    std::unique_ptr<const ALinearOp>&& op,
    std::unique_ptr<const IProj>&& transferOp)
  {
    _operators.push_back(std::move(op));
    _transferOps.push_back(std::move(transferOp));
    _work.push_back(VectorDouble(_operators.back()->getSize()));
    _res.push_back(VectorDouble(_operators.back()->getSize()));
    _rhs_storage.push_back(VectorDouble(_operators.back()->getSize()));
    _err_storage.push_back(VectorDouble(_operators.back()->getSize()));
    _chebys.emplace_back();
    _nLevels++;
  }

  void MultiGridSolver::setCoarseSolver(const MatrixSparse& mat)
  {
    _operators.push_back(std::make_unique<MatrixSparse>(mat));
    _coarseSolver = std::make_unique<const CholeskySparse>(mat);
    _work.push_back(VectorDouble(mat.getNRows()));
    _res.push_back(VectorDouble(mat.getNRows()));
    _rhs_storage.push_back(VectorDouble(mat.getNRows()));
    _err_storage.push_back(VectorDouble(mat.getNRows()));
  }

  void MultiGridSolver::prepare(size_t niterPower)
  {
    _maxEigenValues.clear();
    for (size_t i = 0; i < _nLevels - 1; ++i)
    {
      _maxEigenValues.push_back(
        LH::powerIteration(_operators[i].get(), niterPower));
    }
  }

  VectorDouble MultiGridSolver::vCycle(
    const VectorDouble& rhs,
    const VectorDouble& u,
    size_t lvl) const
  {
    VectorDouble res = u;
    _vCycle(lvl, rhs, res);
    return res;
  }

  void MultiGridSolver::applyLevel(constvect rhs, vect x, size_t lvl) const
  {
    const Chebychev* cheby = &_chebys[lvl];
    double lmax = _maxEigenValues[lvl];
    Id smoothiter = _smoothIter * std::pow(3, lvl);
    cheby->smootherInPlace(
      *_operators[lvl].get(), x, rhs, lmax, _ratiochebmin, _ratiochebmax,
      smoothiter);
  }

  void MultiGridSolver::display() const
  {
    std::cout << "MultiGridSolver with " << _nLevels << " levels." << std::endl;
    for (size_t i = 0; i < _nLevels - 1; ++i)
    {
      std::cout << " Level " << i
                << ": Operator size = " << _operators[i]->getSize()
                << std::endl;
    }
    if (_coarseSolver)
      std::cout << " Coarse solver: " << _coarseSolver->getSize() << " x "
                << _coarseSolver->getSize() << std::endl;
  }

  void MultiGridSolver::_vCycle(size_t lvl, constvect f, vect u) const
  {
    if (lvl == _nLevels - 1)
    {
      _coarseSolver->solve(f, u);
      return;
    }

    // 1. Pre-smoothing
    applyLevel(f, u, lvl);

    // 2. Residual : r = f - A*u
    _operators[lvl]->evalDirect(u, _work[lvl]);

    // subtractInPlace(in1, in2, res) -> res = in2 - in1
    VH::subtractInPlace(_work[lvl], f, _res[lvl]);

    // 3. Restriction
    _transferOps[lvl]->point2mesh(_res[lvl], _rhs_storage[lvl + 1]);

    // 4. Recursion
    _err_storage[lvl + 1].fill(0, _err_storage[lvl + 1].size());
    _vCycle(lvl + 1, _rhs_storage[lvl + 1], _err_storage[lvl + 1]);

    // 5. Prolongation : u = u + P * e
    _transferOps[lvl]->addMesh2point(_err_storage[lvl + 1], u);

    // 6. Post-smoothing
    applyLevel(f, u, lvl);
  }

} // namespace gstlrn
