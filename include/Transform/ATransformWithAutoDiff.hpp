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

#include "Transform/ATransform.hpp"
#ifndef SWIG
#include <boost/math/differentiation/autodiff.hpp>
#endif

// Classe intermédiaire générique utilisant Boost.Autodiff
namespace gstlrn
{
  template<typename Derived>
  class GSTLEARN_EXPORT ATransformWithAutoDiff: public ATransform
  {

  public:
    ATransformWithAutoDiff()
      : ATransform()
    {
    }

    double transform(double x) const override { return evalGeneric<double>(x); }

    double evalJacobian(double h) const override
    {
      // Création variable autodiff pour 1re dérivée
      auto x = boost::math::differentiation::make_fvar<double, 1>(h);

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
