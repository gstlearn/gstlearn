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

#include "LinearOp/ALinearOp.hpp"

#ifndef SWIG

#  include <Eigen/Core>
#  include <Eigen/Dense>
#  include <Eigen/IterativeLinearSolvers>
#  include <Eigen/src/Core/Matrix.h>
#  include <cassert>
// iostream is included here as it is used in Eigen function (std::cerr)
#  include <iostream>
#  include <unsupported/Eigen/IterativeSolvers>

#  define DECLARE_EIGEN_TRAITS(TLinOP)                                                         \
    namespace gstlrn                                                                           \
    {                                                                                          \
    class TLinOP;                                                                              \
    }                                                                                          \
    using Eigen::SparseMatrix;                                                                 \
                                                                                               \
    namespace Eigen                                                                            \
    {                                                                                          \
    namespace internal                                                                         \
    {                                                                                          \
    template<>                                                                                 \
    struct traits<gstlrn::TLinOP>: public Eigen::internal::traits<Eigen::SparseMatrix<double>> \
    {                                                                                          \
    };                                                                                         \
    }                                                                                          \
    }

#  define DECLARE_EIGEN_PRODUCT(TLinOP)                                                                                             \
    template<typename Rhs>                                                                                                          \
    struct Eigen::internal::generic_product_impl<gstlrn::TLinOP, Rhs, Eigen::SparseShape, Eigen::DenseShape, Eigen::GemvProduct>    \
      : Eigen::internal::generic_product_impl_base<gstlrn::TLinOP, Rhs, Eigen::internal::generic_product_impl<gstlrn::TLinOP, Rhs>> \
    {                                                                                                                               \
      typedef typename Product<gstlrn::TLinOP, Rhs>::Scalar Scalar;                                                                 \
      template<typename Dest>                                                                                                       \
      static void scaleAndAddTo(Dest& dst, const gstlrn::TLinOP& lhs, const Rhs& rhs, const Scalar& alpha)                          \
      {                                                                                                                             \
        assert(alpha == Scalar(1) && "scaling is not implemented");                                                                 \
        EIGEN_ONLY_USED_FOR_DEBUG(alpha);                                                                                           \
        lhs.addToDest(rhs, dst);                                                                                                    \
      }                                                                                                                             \
    };
#endif

namespace gstlrn
{
template<typename TLinOP>
class ALinearOpEigenCG: public Eigen::EigenBase<TLinOP>, // No Export because it's a template
                        public virtual ALinearOp         // virtual for ASPDEOp
{
public:
  virtual ~ALinearOpEigenCG() {};

#ifndef SWIG

public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef Id StorageIndex;
  enum
  {
    ColsAtCompileTime    = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor           = false
  };

  Eigen::Index rows() const { return getSize(); }
  Eigen::Index cols() const { return getSize(); }

  template<typename Rhs>
  Eigen::Product<TLinOP, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const
  {
    return Eigen::Product<TLinOP, Rhs, Eigen::AliasFreeProduct>(*(dynamic_cast<const TLinOP*>(this)), x.derived());
  }
#endif
};
} // namespace gstlrn
