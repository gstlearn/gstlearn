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
#include "Basic/Vector.hpp"

class SimuFFTParam;
class DbGrid;
class Model;

class GSTLEARN_EXPORT SimuFFT: public ACalcSimulation
{
public:
  SimuFFT(int nbsimu = 0, int seed = 4324324);
  SimuFFT(const SimuFFT &r) = delete;
  SimuFFT& operator=(const SimuFFT &r) = delete;
  virtual ~SimuFFT();

  int simulate(DbGrid *db,
               Model* model,
               const SimuFFTParam& param,
               int iptr,
               bool verbose = false);
  VectorDouble getChangeSupport(DbGrid *db,
                                Model *model,
                                const SimuFFTParam& param,
                                const VectorDouble& sigma = VectorDouble(),
                                bool verbose = false);

private:
  virtual bool _run() override;

  bool _isValid(Db *db, Model *model);
  void _alloc(DbGrid *db,
              Model *model,
              const SimuFFTParam& param,
              bool verbose = false);
  int _getOptimalEvenNumber(int number, int largeFactor = 11);
  VectorInt _getFactors(int number);
  void _gridDilate(const DbGrid *db,
                   Model *model,
                   const SimuFFTParam& param,
                   bool verbose);
  bool _checkCorrect(Model *model,
                     const VectorVectorDouble& xyz,
                     int ix,
                     int iy,
                     int iz,
                     double percent);
  void _prepar(DbGrid *db,
               Model *model,
               const SimuFFTParam& param,
               bool flag_amplitude,
               bool verbose = false,
               double eps = EPSILON5);
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
  int _ndim;
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
