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
#include "Covariances/ACovFunc.hpp"

class CovContext;
class TurningBandOperate;
class MatrixRectangular;

class GSTLEARN_EXPORT CovMatern : public ACovFunc
{
public:
  CovMatern(const CovContext& ctx);
  CovMatern(const CovMatern &r);
  CovMatern& operator= (const CovMatern &r);
  virtual ~CovMatern();

  virtual String getFormula() const override;
  String         getCovName() const override { return "Matern"; }
  int            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceR() const override { return true; }
  bool           getCompatibleSpaceS() const override { return true; }

  bool   hasParam() const override { return true; }
  double getParMax() const override { return MAX_PARAM; }
  double getScadef() const override;
  bool   hasSpectrumOnSphere() const override { return true; }
  bool   hasSpectrumOnRn() const override { return true; }
  bool   hasMarkovCoeffs() const override { return true; }
  double evaluateSpectrum(double freq) const override;
  void   setMarkovCoeffs(VectorDouble coeffs) override { _markovCoeffs = coeffs;}
  VectorDouble getMarkovCoeffs() const override;
  double getCorrec() const override { return _correc;}
  void   computeCorrec(int ndim) override;
  void   setCorrec(double val) override { _correc = val;}
  void   computeMarkovCoeffs(int dim) override;

  bool isValidForTurningBand() const override { return true; }
  double simulateTurningBand(double t0, TurningBandOperate &operTB) const override;

  bool isValidForSpectral() const override { return true; }
  MatrixRectangular simulateSpectralOmega(int nb) const override;

protected:
  double _evaluateCov(double h) const override;
  VectorDouble _evaluateSpectrumOnSphere(int n, double scale = 1.) const override;

private:
  double _newMatern(double h) const;
  double _oldMatern(double h) const;
private:
  double _correc;
  VectorDouble _markovCoeffs;
};
