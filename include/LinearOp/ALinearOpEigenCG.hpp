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

#include "LinearOp/CGParam.hpp"
#include "LinearOp/LogStats.hpp"
#include "LinearOp/ILinearOpEigenCG.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"

#include "Matrix/VectorEigen.hpp"

#ifndef SWIG
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/src/Core/Matrix.h>
#include <unsupported/Eigen/IterativeSolvers>


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
    lhs.evalDirectEigen(rhs, dst); \
  } \
};
#endif

template<typename TLinOP>
class GSTLEARN_EXPORT ALinearOpEigenCG: public Eigen::EigenBase<TLinOP>,
                                        public ILinearOpEigenCG
{
public:
  ALinearOpEigenCG(const CGParam& params = CGParam());
  ALinearOpEigenCG(const ALinearOpEigenCG &m);
  ALinearOpEigenCG& operator=(const ALinearOpEigenCG &m);
  virtual ~ALinearOpEigenCG();

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

  void evalInverse(const VectorDouble& inv, VectorDouble& outv) const override;
  void evalInverse(const VectorEigen& inv, VectorEigen& outv) const override;

  void evalDirect(const VectorDouble& inv, VectorDouble& outv) const override;
  void evalDirect(const VectorEigen& inv, VectorEigen& outv) const override;

  void setX0(const VectorDouble& x0)  override { _params.setX0(x0); }
  void mustShowStats(bool status)     override { _logStats.mustShowStats(status); }

  const LogStats& getLogStats() const override { return _logStats; }

#ifndef SWIG
  virtual void evalInverseEigen(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const;
  virtual void evalDirectEigen(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const;
#endif

protected:
  virtual void _evalDirectEigen(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const = 0;

private:
  double _prod(const VectorDouble& x, const VectorDouble& y) const;

private:
  CGParam _params;

protected:
  LogStats _logStats;
};

#ifndef SWIG

template <typename TLinOP>
ALinearOpEigenCG<TLinOP>::ALinearOpEigenCG(const CGParam& params)
  : _params(params),
    _logStats()
{
}

template <typename TLinOP>
ALinearOpEigenCG<TLinOP>::ALinearOpEigenCG(const ALinearOpEigenCG &m)
  : _params(m._params),
    _logStats(m._logStats)
{
}

template <typename TLinOP>
ALinearOpEigenCG<TLinOP>& ALinearOpEigenCG<TLinOP>::operator=(const ALinearOpEigenCG<TLinOP> &m)
{
  if (this != &m)
  {
    _params = m._params;
    _logStats = m._logStats;
  }
  return *this;
}

template <typename TLinOP>
ALinearOpEigenCG<TLinOP>::~ALinearOpEigenCG() 
{
  _logStats.statsShow();
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = Q * 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
template <typename TLinOP>
void ALinearOpEigenCG<TLinOP>::evalDirect(const VectorDouble& inv, VectorDouble& outv) const
{
  try
  {
    Eigen::Map<const Eigen::VectorXd> myInv(inv.data(), inv.size());
    Eigen::VectorXd myOut;
    // Assume outv has the good size
    _evalDirectEigen(myInv, myOut);
    Eigen::Map<Eigen::VectorXd>(outv.data(), outv.size()) = myOut;
  }
  catch(const std::string& str)
  {
    // TODO : Check if std::exception can be used
    messerr("%s", str.c_str());
  }
}


/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = Q * 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
template<typename TLinOP>
void ALinearOpEigenCG<TLinOP>::evalDirect(const VectorEigen& inv,
                                          VectorEigen& outv) const
{
  evalDirectEigen(inv.getVector(), outv.getVector());
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = Q * 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
template<typename TLinOP>
void ALinearOpEigenCG<TLinOP>::evalDirectEigen(const Eigen::VectorXd& inv,
                                               Eigen::VectorXd& outv) const
{
  try
  {
    _evalDirectEigen(inv, outv);
  }
  catch (const std::string& str)
  {
    // TODO : Check if std::exception can be used
    messerr("%s", str.c_str());
  }
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = Q^{-1} * 'inv' by Conjugate Gradient
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
template <typename TLinOP>
void ALinearOpEigenCG<TLinOP>::evalInverse(const VectorDouble &inv, VectorDouble &outv) const
{
  Eigen::Map<const Eigen::VectorXd> myInv(inv.data(), inv.size());
  Eigen::VectorXd myOut;
  // Assume outv has the good size
  evalInverseEigen(myInv, myOut);
  Eigen::Map<Eigen::VectorXd>(outv.data(), outv.size()) = myOut;
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = Q^{-1} * 'inv' by Conjugate Gradient
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
template<typename TLinOP>
void ALinearOpEigenCG<TLinOP>::evalInverse(const VectorEigen& inv,
                                           VectorEigen& outv) const
{
  evalInverseEigen(inv.getVector(), outv.getVector());
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = Q^{-1} * 'inv' by Conjugate Gradient
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
template<typename TLinOP>
void ALinearOpEigenCG<TLinOP>::evalInverseEigen(const Eigen::VectorXd& inv,
                                                Eigen::VectorXd& outv) const
{
  Eigen::ConjugateGradient<TLinOP,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::IdentityPreconditioner> cg;
  cg.compute(*this);
  outv = cg.solve(inv);
}

/*****************************************************************************/
/*!
**  Returns the scalar product between 'x' and 'y'
**
** \param[in]  x      First array
** \param[in]  y      Second array
**
*****************************************************************************/
template <typename TLinOP>
double ALinearOpEigenCG<TLinOP>::_prod(const VectorDouble &x, const VectorDouble &y) const
{
  double prod = 0.;
  for (int i=0, n = getSize(); i<n; i++)
    prod += x[i] * y[i];
  return prod;
}

#endif