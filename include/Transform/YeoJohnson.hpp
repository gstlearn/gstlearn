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

#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Transform/ATransform.hpp"
#include "Transform/ATransformWithAutoDiff.hpp"

#include "gstlearn_export.hpp"

namespace gstlrn
{

  class GSTLEARN_EXPORT YeoJohnson:
#ifndef SWIG
    public ATransformWithAutoDiff<YeoJohnson>
#else
    public ATransform
#endif
  {
  public:
    YeoJohnson(double lambda);
    YeoJohnson(const YeoJohnson& r);
    YeoJohnson& operator=(const YeoJohnson& r);
    virtual ~YeoJohnson() = default;
    IMPLEMENT_CLONING(YeoJohnson)

    bool hasParameters() const override { return true; }

    void _printParams(std::stringstream& sstr,
                      const AStringFormat* strfmt) const override;
    double inverseTransform(double x) const override;
    VectorDouble getParams() const override;

    double getLambdaValue() const { return _lambda.getValue(); }

    void setLambdaValue(double lambda) { _lambda.setValue(lambda); }
#ifndef SWIG
    void initParams(double min = 0., double max = INF) override;
#endif
    void appendParams(ListParams& listParams) override;
    static YeoJohnson* create(double lambda);

    String getName() const override { return "Yeo-Johnson"; }

    void setK(double K) { _K = K; }

    double getSaturation() const { return _K; }

    template<typename T>
    T evalImpl(T h) const
    {
      const double lambda = _lambda.getValue();
      const double eps = 1e-6; // marge de sécurité

      if (h >= 0)
      {
        if (std::abs(lambda) < 1e-12) return exp(h) - 1.0;

        T base = 1.0 + (lambda * h);
        T y;

        if (lambda > 0)
        {
          // Cas standard
          y = pow(base, 1.0 / lambda) - 1.0;
        }
        else
        {
          // λ < 0 : risque d’explosion quand base → 0+
          if (base > eps)
          {
            y = pow(base, 1.0 / lambda) - 1.0;
          }
          else
          {
            // Zone critique : on fige la valeur à eps
            y = pow(eps, 1.0 / lambda) - 1.0;
          }
        }

        // Saturation lisse pour éviter les valeurs extrêmes
        return (_K * tanh(y / _K));
      }

      // Cas h < 0
      if (std::abs(lambda - 2.0) < 1e-12) return 1.0 - exp(-h);

      T base_neg = 1.0 - ((2.0 - lambda) * h);
      T y = 1.0 - pow(base_neg, 1.0 / (2.0 - lambda));

      // Saturation lisse également côté négatif
      return (_K * tanh(y / _K));
    }

  private:
    ParamInfo _lambda;
    double _K = 1e6; // Valeur de saturation pour la zone critique
  };

} // namespace gstlrn
