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
#include <Eigen/IterativeLinearSolvers>
#include "LinearOp/APreconditioner.hpp"

namespace gstlrn
{
  
class IdentityPreconditioner : public APreconditioner
{
  public:

    IdentityPreconditioner() {}
    template<typename MatrixType>
    explicit IdentityPreconditioner(const MatrixType& ) {}

    template<typename MatrixType>
    IdentityPreconditioner& analyzePattern(const MatrixType& ) { return *this; }

    template<typename MatrixType>
    IdentityPreconditioner& factorize(const MatrixType& ) { return *this; }

    template<typename MatrixType>
    IdentityPreconditioner& compute(const MatrixType& ) { return *this; }

    template<typename Rhs>
    inline const Rhs& solve(const Rhs& b) const { return b;  }

    Eigen::ComputationInfo info() { return Success; }

    private:
    mutable Eigen::VectorXd _c;
    Eigen::ComputationInfo Success = Eigen::Success;
};


}   