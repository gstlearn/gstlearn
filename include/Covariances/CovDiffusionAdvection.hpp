/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Arrays/Array.hpp"
#include "Covariances/CovContext.hpp"

#include <complex>

class CovAniso;
class GSTLEARN_EXPORT CovDiffusionAdvection
{

public:
  CovDiffusionAdvection();
  CovDiffusionAdvection(const CovDiffusionAdvection& r);
  CovDiffusionAdvection& operator=(const CovDiffusionAdvection& r);
  virtual ~CovDiffusionAdvection();

  static CovDiffusionAdvection* create(CovAniso* markovL = nullptr,
                                       CovAniso* markovR = nullptr,
                                       double scaleTime = 1.,
                                       VectorDouble vel = VectorDouble(2),
                                       double sigma2 = 1.);

  const CovAniso* getMarkovL() const { return _markovL; }
  void setMarkovL(const CovAniso *markovL) { _markovL = markovL; }
  const CovAniso* getMarkovR() const { return _markovR; }
  void setMarkovR(const CovAniso *markovR) { _markovR = markovR; }
  double getScaleTime() const { return _scaleTime; }
  void setScaleTime(double scaleTime) { _scaleTime = scaleTime; }
  double getSigma2() const { return _sigma2; }
  void setSigma2(double sigma2) { _sigma2 = sigma2; }
  const VectorDouble& getVel() const { return _vel; }
  void setVel(const VectorDouble &vel) { _vel = vel; }
  double getGlobalCorrec() const { return _globalCorrec; }
  const CovAniso* getSpatialTrace() const { return _spatialTrace; }

  bool isNoneMarkovL() const {return _markovL==nullptr;}
  bool isNoneMarkovR() const {return _markovR==nullptr;}

  std::complex<double> evalSpatialSpectrum(VectorDouble freq, double time) const;
  Array evalCovFFT(const VectorDouble& hmax,double time = 0, int N = 128) const;

private:

  void _init();
  void _computeSpatialTrace();

  const CovAniso*    _markovL;
  const CovAniso*    _markovR;
  double             _scaleTime;
  VectorDouble       _vel;
  double             _sigma2;
  double             _globalCorrec;
  CovAniso*          _spatialTrace;
  CovContext         _ctxt;
  bool               _destroyMarkovL;
  bool               _destroyMarkovR;
  bool               _markovRdefined;
  bool               _markovLdefined;
};
