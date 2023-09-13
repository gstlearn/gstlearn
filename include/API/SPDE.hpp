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

#include "Enum/ESPDECalcMode.hpp"

#include "Basic/NamingConvention.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"

class ShiftOpCs;
class Db;
class DbGrid;
class PrecisionOp;
class PrecisionOpCs;
class Model;
class MeshETurbo;

/**
 * The SPDE class provides the SPDE implementation of a univariate model defined by
 * the sum of a nugget effect and Matern's models. Its main objectives are:
 * - point kriging with or without linear drifts (SK, OK, KED, UK)
 * - point simulations, conditional or non conditional
 * - evaluation of the log likelihood of the model, in order to estimate the parameter using maximum likelihood
 */
class GSTLEARN_EXPORT SPDE
{
public:
  SPDE(Model *model,
       const Db* domain,
       const Db* data = nullptr,
       const ESPDECalcMode& calcul = ESPDECalcMode::fromKey("SIMUCOND"),
       const AMesh* mesh = nullptr,
       bool verbose = false);
  SPDE(const SPDE& r) = delete;
  SPDE& operator=(const SPDE& r) = delete;
  virtual ~SPDE();

  static SPDE* create(Model *model,
                      const Db *domain,
                      const Db *data = nullptr,
                      const ESPDECalcMode& calcul = ESPDECalcMode::fromKey("SIMUCOND"),
                      const AMesh* mesh = nullptr,
                      bool verbose = false);

  int compute(Db *dbout,
              int nbsimu = 1,
              int seed = 131351,
              const NamingConvention &namconv = NamingConvention("spde"));

  double computeLogDet(int nbsimu = 1,int seed = 1234) const;
  double computeQuad() const;
  double computeLogLike(int nbsimu = 1, int seed = 131323) const;
  double computeProfiledLogLike(int nbsimu = 1, int seed = 131323) const;
  VectorDouble getCoeffs();

  void setDriftCoeffs(VectorDouble coeffs);
  void setEps(double eps) { _eps = eps; }
  void setNIterMax(int nitermax) { _nIterMax = nitermax; }

  const PrecisionOpCs* getPrecisionOp(int i = 0) const  { return (PrecisionOpCs*)_pilePrecisions[i];}
  const ProjMatrix* getProjMatrix(int i = 0) const  { return _pileProjMatrix[i];}
  const PrecisionOpMultiConditional* getPrecisionKriging() const { return _precisionsKrig;}
  const AMesh* getKrigingMeshing(int i = 0) const { return _meshingKrig[i];}
  const AMesh* getSimuMeshing(int i = 0) const { return _meshingSimu[i]; }
  const Db* getData() const {return  _data;}

  void setRefineK(int refineK)        { _refineK = refineK; }
  void setRefineS(int refineS)        { _refineS = refineS; }
  void setBorder(int border)          { _border = border; }
  void setEpsNugget(double epsNugget) { _epsNugget = epsNugget; }

private:
  int _init(Model *model,
            const Db *domain,
            const Db *data = nullptr,
            const ESPDECalcMode &calcul = ESPDECalcMode::fromKey("SIMUCOND"),
            const AMesh *mesh = nullptr,
            bool verbose = false);
  void _centerByDrift(const VectorDouble& dataVect,int ivar=0,bool useSel=true) const;
  void _computeDriftCoeffs() const;
  void _purge();
  bool _performSimulation() const;
  bool _performKriging() const;
  void _computeLk() const;
  void _computeKriging() const;
  void _computeSimuNonCond() const;
  void _computeSimuCond() const;
  void _addNuggetOnResult(VectorDouble &result);
  void _addDrift(Db* db, VectorDouble &result, int ivar = 0, bool useSel = true);
  bool _useCholesky() const;

private:
  const Db*_data;
  ESPDECalcMode _calcul;
  int _refineK;
  int _refineS;
  int _border;
  PrecisionOpMultiConditional* _precisionsKrig;
  PrecisionOpMultiConditional* _precisionsSimu;
  std::vector<PrecisionOp*>    _pilePrecisions; // Dimension: number of valid covariances
  std::vector<ProjMatrix*>     _pileProjMatrix; // Dimension: number of valid covariances
  std::vector<const AMesh*>    _meshingSimu;    // Dimension: number of valid covariances
  std::vector<const AMesh*>    _meshingKrig;    // Dimension: number of valid covariances
  mutable VectorDouble         _driftCoeffs;
  Model*                       _model;
  mutable VectorVectorDouble   _workingKrig;     // Number of Mesh apices * Number of valid covariances
  mutable VectorVectorDouble   _workingSimu;     // Number of Mesh apices * Number of valid covariances
  mutable VectorDouble         _workingData;     // Number of valid data
  mutable VectorDouble         _workingDataInit; // Number of valid data
  std::vector<ProjMatrix*>     _projOnDbOut;
  VectorInt                    _adressesICov;
  double _nugget;
  VectorVectorDouble _driftTab;
  bool _requireCoeffs;
  mutable bool _isCoeffsComputed;
  bool _deleteMesh;

  // Parameters specific invertion using Conjugate Gradient (used for Kriging)
  int _nIterMax;
  double _eps;
  double _epsNugget; // Additional amount of nugget specified as the percentage of total sill (nugget excluded)
};

GSTLEARN_EXPORT int krigingSPDE(Db *dbin,
                                Db *dbout,
                                Model *model,
                                bool flag_est = true,
                                bool flag_std = false,
                                bool flag_varz = false,
                                const AMesh* mesh = nullptr,
                                int refineK = 11,
                                int border = 8,
                                double epsNugget = 1.e-2,
                                bool verbose = false,
                                const NamingConvention &namconv = NamingConvention("KrigingSPDE"));
GSTLEARN_EXPORT int simulateSPDE(Db *dbin,
                                 Db *dbout,
                                 Model *model,
                                 int nbsimu = 1,
                                 const AMesh *mesh = nullptr,
                                 int refineK = 11,
                                 int refineS = 18,
                                 int border = 8,
                                 int seed = 121423,
                                 double epsNugget = 1.e-2,
                                 bool verbose = false,
                                 const NamingConvention &namconv = NamingConvention("SimuSPDE"));
