/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "geoslib_define.h"

#include "Simulation/ACalcSimulation.hpp"
#include "Basic/Plane.hpp"
#include "Basic/VectorNumT.hpp"

class SimuFFTParam;
class DbGrid;
class Model;

class GSTLEARN_EXPORT CalcSimuFFT: public ACalcSimulation
{
public:
  CalcSimuFFT(int nbsimu = 0, bool verbose = false, int seed = 4324324);
  CalcSimuFFT(const CalcSimuFFT &r) = delete;
  CalcSimuFFT& operator=(const CalcSimuFFT &r) = delete;
  virtual ~CalcSimuFFT();

  void setParam(const SimuFFTParam &param) { _param = param; }
  void setVerbose(bool verbose) { _verbose = verbose; }
  VectorDouble changeSupport(const VectorDouble &sigma);

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  bool _simulate();
  void _alloc();
  int _getOptimalEvenNumber(int number, int largeFactor = 11);
  VectorInt _getFactors(int number);
  void _gridDilate();
  bool _checkCorrect(const VectorVectorDouble& xyz,
                     int ix,
                     int iy,
                     int iz,
                     double percent);
  void _prepar(bool flag_amplitude, double eps = EPSILON5);
  void _defineRandom();
  void _setVariance(int ix, int iy, int iz);
  void _defineSymmetry(void);
  void _defineSym1();
  void _defineSym2(int iz0);
  void _defineSym3();
  void _setZero(int ix, int iy, int iz);
  void _setConjugate(int ix, int iy, int iz, int jx, int jy, int jz);
  void _final(DbGrid *db, int iad);
  double _support(double sigma);
  double _support1(double sigma);
  double _support2(double sigma);
  double _support3(double sigma);
  double _rhoSigma(double sigma, int ix, int iy, int iz);

private:
  int _iattOut;
  bool _verbose;
  SimuFFTParam _param;
  int _nxyz;
  VectorInt _nx;
  VectorInt _shift;
  VectorInt _dims;
  VectorInt _dim2;
  int _sizes_alloc;
  VectorDouble _cmat;
  VectorDouble _rnd;
  VectorDouble _u;
  VectorDouble _v;
};

GSTLEARN_EXPORT int simfft(DbGrid *db,
                           Model *model,
                           SimuFFTParam& param,
                           int nbsimu = 1,
                           int seed = 432431,
                           int verbose = false,
                           const NamingConvention& namconv = NamingConvention("FFT"));
GSTLEARN_EXPORT VectorDouble getChangeSupport(DbGrid *db,
                                              Model *model,
                                              const SimuFFTParam &param,
                                              const VectorDouble &sigma = VectorDouble(),
                                              int seed = 14333,
                                              bool verbose = false);
