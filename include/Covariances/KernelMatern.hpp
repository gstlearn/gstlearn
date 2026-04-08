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
#include "Enum/ESimuType.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{

  class CovContext;
  class TurningBandOperate;
  class MatrixDense;

  class GSTLEARN_EXPORT KernelMatern: public AKernel
  {
  public:
    KernelMatern(const CovContext& ctx);
    KernelMatern(const KernelMatern& r);
    KernelMatern& operator=(const KernelMatern& r);
    virtual ~KernelMatern();

    String getFormula() const override;

    String getCovName() const override { return "Matern"; }

    Id getMinOrder() const override { return -1; }

    bool getCompatibleSpaceR() const override { return true; }

    bool getCompatibleSpaceS() const override { return true; }

    bool hasParam() const override { return true; }

    double getParMax() const override { return MAX_PARAM; }

    double getScadef() const override;

    bool hasMarkovCoeffs() const override { return true; }

    void setMarkovCoeffs(const VectorDouble& coeffs) override
    {
      _markovCoeffs = coeffs;
    }

    VectorDouble getMarkovCoeffs() const override;
    void computeMarkovCoeffs(Id dim) override;

    double getCorrec() const override { return _correc; }

    void computeCorrec(Id ndim) override;

    void setCorrec(double val) override { _correc = val; }

    double
      simulateTurningBand(double t0, TurningBandOperate& operTB) const override;

    bool isValidForSimulation(const ESimuType& simuType) const override
    {
      return (simuType == ESimuType::SPECTRAL || simuType == ESimuType::TB);
    }

    double evaluateSpectrum(double freq) const override;
    MatrixDense simulateSpectralOmega(Id nb) const override;

    bool needCorrec() const override { return true; }

  protected:
    double _evaluateCov(double h) const override;
    double _evaluateCovGeneric(double h) const;

    double _evaluateCovFirstDerivative(double h) const override;
    VectorDouble
      _evaluateSpectrumOnSphere(Id n, double scale = 1., bool flagScale = true)
        const override;
    void _setParam(double param, Id ipar = 0) override;

  private:
    using MaternFunc = double (*)(double, Id);
    static double _besselK(double nu, double h);
    static double _evalExp(double h, Id p = 0);
    static double _evalNu15(double h, Id p = 1);
    static double _evalNu25(double h, Id p = 2);
    static double _evaluateCovIntegerPlusOneHalf(double h, Id p);

  private:
    MaternFunc _maternFunc = nullptr;
    double _correc;
    VectorDouble _markovCoeffs;
  };

} // namespace gstlrn
