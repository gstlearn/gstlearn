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
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
DISABLE_WARNING_POP
#endif

class css; /// TODO : Dependency to csparse to be removed
class csn;
class MatrixSparse;

class CholeskyEigenCG;

// Eigen stuff for conjugate gradient
/////////////////////////////////////////////////////////////////////////////////
namespace Eigen {
namespace internal {
  // CholeskyEigenCG looks-like a SparseMatrix, so let's inherits its traits:
  template<>
  struct traits<CholeskyEigenCG> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
  {};
}
}
/////////////////////////////////////////////////////////////////////////////////

class GSTLEARN_EXPORT CholeskyEigenCG: public ALinearOp, public Eigen::EigenBase<CholeskyEigenCG> 
{
public:
  CholeskyEigenCG(const MatrixSparse* mat);
  CholeskyEigenCG(const CholeskyEigenCG &m) = delete;
  CholeskyEigenCG& operator=(const CholeskyEigenCG &m) = delete;
  virtual ~CholeskyEigenCG();

  // Eigen stuff for conjugate gradient
  /////////////////////////////////////////////////////////////////////////////////
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };
  Index rows() const { return getSize(); }
  Index cols() const { return getSize(); }

  template<typename Rhs>
  Eigen::Product<CholeskyEigenCG,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<CholeskyEigenCG,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }
  /////////////////////////////////////////////////////////////////////////////////

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
  css *_S; // CholeskyEigenCG decomposition (for Old-style Csparse storage)
  csn *_N; // CholeskyEigenCG decomposition (for Old-style Csparse storage)
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > _cholSolver; // (for Eigen library storage)
#endif
  const MatrixSparse* _matCS; // Stored by compliance with ALinearOp. Not to be deleted
};


// Implementation of ALinearOpEigenCG * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {
 
  template<typename Rhs>
  struct generic_product_impl<CholeskyEigenCG, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<CholeskyEigenCG,Rhs,generic_product_impl<CholeskyEigenCG,Rhs> >
  {
    typedef typename Product<CholeskyEigenCG,Rhs>::Scalar Scalar;
 
    template<typename Dest>
    static void scaleAndAddTo(Dest& dst, const CholeskyEigenCG& lhs, const Rhs& rhs, const Scalar& alpha)
    {
      // This method should implement "dst += alpha * lhs * rhs" inplace,
      // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
      assert(alpha==Scalar(1) && "scaling is not implemented");
      EIGEN_ONLY_USED_FOR_DEBUG(alpha);
 
      // Matrix free product used by Eigen Conjugate Gradient
      // TODO : Too much memory copies
      VectorDouble inv(rhs.data(), rhs.data() + rhs.size());
      VectorDouble outv(inv.size());
      lhs.evalDirect(inv, outv);
      for(Index i=0; i<lhs.cols(); ++i)
        dst[i] = outv[i];
    }
  };
 
}
}
