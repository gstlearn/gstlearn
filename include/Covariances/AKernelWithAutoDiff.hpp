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
#include "Covariances/AKernel.hpp"
#ifndef SWIG
#  include <boost/math/differentiation/autodiff.hpp>
#endif

// Classe intermédiaire générique utilisant Boost.Autodiff
namespace gstlrn
{
template<typename Derived>
class AKernelWithAutoDiff: public AKernel
{
public:
  AKernelWithAutoDiff(const ECov& type, const CovContext& ctxt)
    : AKernel(type, ctxt)
  {
  }
  AKernelWithAutoDiff(const AKernelWithAutoDiff& r)
    : AKernel(r)
  {
  }
  AKernelWithAutoDiff& operator=(const AKernelWithAutoDiff& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }
  double _evaluateCov(double h) const override
  {
    return evalGeneric<double>(h); // Appel version double
  }

  double _evaluateCovFirstDerivative(double h) const override
  {
    // Création variable autodiff pour 1re dérivée
    boost::math::differentiation::autodiff_fvar<double, 1> x = boost::math::differentiation::make_fvar<double, 1>(h);

    // Appel version générique (autodiff)
    auto y = evalGeneric(x);

    return y.derivative(1); // ∂f/∂x
  }

protected:
  // Implémentation spécifique dans la classe dérivée
  template<typename T>
  T evalGeneric(T h) const
  {
    // Redirige vers l'implémentation de la classe fille
    return static_cast<const Derived*>(this)->evalImpl(h);
  }
};
} // namespace gstlrn