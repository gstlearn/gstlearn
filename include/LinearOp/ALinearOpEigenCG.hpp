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
#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"

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
    lhs.evalDirect(rhs, dst); \
  } \
};

template <typename TLinOP>
class GSTLEARN_EXPORT ALinearOpEigenCG : public Eigen::EigenBase<TLinOP> {

public:
  ALinearOpEigenCG(const CGParam params = CGParam());
  ALinearOpEigenCG(const ALinearOpEigenCG &m);
  ALinearOpEigenCG& operator=(const ALinearOpEigenCG &m);
  virtual ~ALinearOpEigenCG();

public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };

  Eigen::Index rows() const { return getSize(); }
  Eigen::Index cols() const { return getSize(); }
 
  template<typename Rhs>
  Eigen::Product<TLinOP,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<TLinOP,Rhs,Eigen::AliasFreeProduct>(*(dynamic_cast<const TLinOP*>(this)), x.derived());
  }

  virtual void evalInverse(const VectorDouble& inv, VectorDouble& outv) const;
  virtual void evalInverse(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const;
  virtual int getSize() const = 0;

  void evalDirect(const VectorDouble& inv, VectorDouble& outv) const;
  void evalDirect(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const;

  void setX0(const VectorDouble& x0) { _params.setX0(x0); }
  void mustShowStats(bool status) { _logStats.mustShowStats(status); }

  const LogStats& getLogStats() const { return _logStats; }

protected:
  virtual void _evalDirect(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const = 0;

private:
  double _prod(const VectorDouble& x, const VectorDouble& y) const;

private:
  CGParam _params;

protected:
  LogStats _logStats;
};

template <typename TLinOP>
ALinearOpEigenCG<TLinOP>::ALinearOpEigenCG(const CGParam params)
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
    Eigen::Map<Eigen::VectorXd> myOut(outv.data(), outv.size());
    _evalDirect(myInv,myOut);
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
template <typename TLinOP>
void ALinearOpEigenCG<TLinOP>::evalDirect(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const
{
  try
  {
    _evalDirect(inv,outv);
  }
  catch(const std::string& str)
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
  //Eigen::Map<Eigen::VectorXd> myOut(outv.data(), outv.size());
  Eigen::VectorXd myOut;
  evalInverse(myInv, myOut);
  outv.assign(myOut.data(), myOut.data()+myOut.size());
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
void ALinearOpEigenCG<TLinOP>::evalInverse(const Eigen::VectorXd &inv, Eigen::VectorXd &outv) const
{
  Eigen::ConjugateGradient<TLinOP, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg;
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