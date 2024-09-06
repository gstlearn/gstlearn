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
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/src/Core/Matrix.h>
#include <unsupported/Eigen/IterativeSolvers>
#include <cassert>


#define DECLARE_EIGEN_TRAITS(TLinOP) \
class TLinOP; \
using Eigen::SparseMatrix; \
 \
namespace Eigen { \
namespace internal { \
  template<> \
  struct traits<TLinOP> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> > \
  {}; \
} \
}

#define DECLARE_EIGEN_PRODUCT(TLinOP) \
template<typename Rhs> \
struct Eigen::internal::generic_product_impl<TLinOP, Rhs, Eigen::SparseShape, Eigen::DenseShape, Eigen::GemvProduct> \
: Eigen::internal::generic_product_impl_base<TLinOP, Rhs, Eigen::internal::generic_product_impl<TLinOP,Rhs> > \
{ \
  typedef typename Product<TLinOP,Rhs>::Scalar Scalar; \
  template<typename Dest> \
  static void scaleAndAddTo(Dest& dst, const TLinOP& lhs, const Rhs& rhs, const Scalar& alpha) \
  { \
    assert(alpha==Scalar(1) && "scaling is not implemented"); \
    EIGEN_ONLY_USED_FOR_DEBUG(alpha); \
    lhs.addToDest(rhs, dst); \
  } \
};
#endif

template<typename TLinOP>
class ALinearOpEigenCG : public Eigen::EigenBase<TLinOP>, // No Export because it's a template
                         public ALinearOp
{
public:
  virtual ~ALinearOpEigenCG() {};

#ifndef SWIG
public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum
  {
    ColsAtCompileTime    = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor           = false
  };
  
  Eigen::Index rows() const { return getSize(); }
  Eigen::Index cols() const { return getSize(); }
 
  template<typename Rhs>
  Eigen::Product<TLinOP,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<TLinOP,Rhs,Eigen::AliasFreeProduct>(*(dynamic_cast<const TLinOP*>(this)), x.derived());
  }
#endif
};

