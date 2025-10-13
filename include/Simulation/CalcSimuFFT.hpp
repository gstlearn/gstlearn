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

#include "geoslib_define.h"

#include "Basic/Plane.hpp"
#include "Basic/VectorNumT.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/SimuFFTParam.hpp"

namespace gstlrn
{

class SimuFFTParam;
class DbGrid;
class Model;

class GSTLEARN_EXPORT CalcSimuFFT: public ACalcSimulation
{
public:
  CalcSimuFFT(Id nbsimu = 0, bool verbose = false, Id seed = 4324324);
  CalcSimuFFT(const CalcSimuFFT& r)            = delete;
  CalcSimuFFT& operator=(const CalcSimuFFT& r) = delete;
  virtual ~CalcSimuFFT();

  void setParam(const SimuFFTParam& param) { _param = param; }
  void setVerbose(bool verbose) { _verbose = verbose; }
  VectorDouble changeSupport(const VectorDouble& sigma);

private:
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;
  bool _postprocess() override;
  void _rollback() override;

  bool _simulate();
  void _alloc();
  static Id _getNOptimalEven(Id number, Id largeFactor = 11);
  static VectorInt _getFactors(Id number);
  void _gridDilate();
  bool _checkCorrect(const VectorVectorDouble& xyz,
                     Id ix,
                     Id iy,
                     Id iz,
                     double percent);
  void _prepar(bool flag_amplitude, double eps = EPSILON5);
  void _defineRandom();
  void _setVariance(Id ix, Id iy, Id iz);
  void _defineSymmetry(void);
  void _defineSym1();
  void _defineSym2(Id iz0);
  void _defineSym3();
  void _setZero(Id ix, Id iy, Id iz);
  void _setConjugate(Id ix, Id iy, Id iz, Id jx, Id jy, Id jz);
  void _final(DbGrid* db, Id iad);
  double _support(double sigma);
  double _support1(double sigma);
  double _support2(double sigma);
  double _support3(double sigma);
  double _rhoSigma(double sigma, Id ix, Id iy, Id iz);

private:
  Id _iattOut;
  bool _verbose;
  SimuFFTParam _param;
  Id _nxyz;
  VectorInt _nx;
  VectorInt _shift;
  VectorInt _dims;
  VectorInt _dim2;
  Id _sizes_alloc;
  VectorDouble _cmat;
  VectorDouble _rnd;
  VectorDouble _u;
  VectorDouble _v;
};

GSTLEARN_EXPORT Id simfft(DbGrid* db,
                           ModelGeneric* model,
                           SimuFFTParam& param,
                           Id nbsimu                      = 1,
                           Id seed                        = 432431,
                           Id verbose                     = false,
                           const NamingConvention& namconv = NamingConvention("FFT"));
GSTLEARN_EXPORT VectorDouble getChangeSupport(DbGrid* db,
                                              ModelGeneric* model,
                                              const SimuFFTParam& param,
                                              const VectorDouble& sigma = VectorDouble(),
                                              Id seed                  = 14333,
                                              bool verbose              = false);
} // namespace gstlrn