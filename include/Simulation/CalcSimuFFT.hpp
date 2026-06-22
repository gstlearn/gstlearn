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

#include "Basic/VectorNumT.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Simulation/ACalcSimuGaussian.hpp"
#include "Simulation/SimuFFTParam.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{

  class SimuFFTParam;
  class DbGrid;
  class Model;

  class GSTLEARN_EXPORT CalcSimuFFT: public ACalcSimuGaussian
  {
  public:
    CalcSimuFFT(Id nbsimu = 0, bool verbose = false, Id seed = 4324324);
    CalcSimuFFT(const CalcSimuFFT& r) = delete;
    CalcSimuFFT& operator=(const CalcSimuFFT& r) = delete;
    virtual ~CalcSimuFFT();

    void setParam(const SimuFFTParam& param) { _param = param; }

    VectorDouble changeSupport(const VectorDouble& sigma);

  private:
    bool _check() override;
    bool _preprocess() override;
    void _rollback() override;

    bool _initializeSimulations() override;
    bool _simulate(Id isimu) override;
    void
      _compute(Db* db, const VectorBool& activeArray, VectorVectorDouble& tab)
        override;

    void _allocate();
    static Id _getNOptimalEven(Id number, Id largeFactor = 11);
    static VectorInt _getFactors(Id number);
    void _gridDilate();
    bool _checkCorrect(
      const VectorVectorDouble& xyz,
      Id ix,
      Id iy,
      Id iz,
      double percent);
    void _prepar(bool flag_amplitude, double eps = EPSILON5);
    void _setVariance(Id ix, Id iy, Id iz);
    void _defineSymmetry(void);
    void _defineSym1();
    void _defineSym2(Id iz0);
    void _defineSym3();
    void _setZero(Id ix, Id iy, Id iz);
    void _setConjugate(Id ix, Id iy, Id iz, Id jx, Id jy, Id jz);

    double _support(double sigma);
    double _support1(double sigma);
    double _support2(double sigma);
    double _support3(double sigma);
    double _rhoSigma(double sigma, Id ix, Id iy, Id iz);

  private:
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
    NeighUnique* _neigh; // This is a useless member, but needed for _check()
  };

} // namespace gstlrn
