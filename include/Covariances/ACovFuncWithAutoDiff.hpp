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
#include "Covariances/ACovFunc.hpp"
#ifndef SWIG
#  include <boost/math/differentiation/autodiff.hpp>
#endif

// Classe intermédiaire générique utilisant Boost.Autodiff
namespace gstlrn
{
template<typename Derived>
class ACovFuncWithAutoDiff: public ACovFunc
{
public:
  ACovFuncWithAutoDiff(const ECov& type, const CovContext& ctxt)
    : ACovFunc(type, ctxt)
  {
  }
  ACovFuncWithAutoDiff(const ACovFuncWithAutoDiff& r)
    : ACovFunc(r)
  {
  }
  ACovFuncWithAutoDiff& operator=(const ACovFuncWithAutoDiff& r)
  {
    if (this != &r)
    {
      ACovFunc::operator=(r);
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