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
#include "geoslib_define.h"

#include "Enum/ESPDECalcMode.hpp"
#include "Basic/NamingConvention.hpp"
#include "API/SPDEParam.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include <vector>

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
       int useCholesky = -1,
       const SPDEParam& params = SPDEParam(),
       bool verbose = false,
       bool showStats = false);
  SPDE(const SPDE& r) = delete;
  SPDE& operator=(const SPDE& r) = delete;
  virtual ~SPDE();

  static SPDE* create(Model *model,
                      const Db *domain,
                      const Db *data = nullptr,
                      const ESPDECalcMode& calcul = ESPDECalcMode::fromKey("SIMUCOND"),
                      const AMesh* mesh = nullptr,
                      int useCholesky = -1,
                      const SPDEParam& params = SPDEParam(),
                      bool verbose = false,
                      bool showStats = false);

  int compute(Db *dbout,
              int nbsimu = 1,
              const NamingConvention &namconv = NamingConvention("spde"));

  double computeLogDet(int nbsimu = 1) const;
  double computeQuad() const;
  double computeLogLikelihood(int nbsimu = 1) const;
  VectorDouble getCoeffs();

  void setDriftCoeffs(const VectorDouble& coeffs);

  const PrecisionOpCs* getPrecisionOpCs(int i = 0) const  { return (PrecisionOpCs*) _pilePrecisions[i];}
  const ProjMatrix* getProjMatrix(int i = 0) const  { return _pileProjMatrix[i];}
  const PrecisionOpMultiConditional* getPrecisionKrig() const { return _precisionsKrig;}
  const AMesh* getMeshingKrig(int i = 0) const { return _meshingKrig[i];}
  const AMesh* getMeshingSimu(int i = 0) const { return _meshingSimu[i]; }
  const Db* getData() const {return  _data;}

private:
  int _init(const Db *domain,
            const AMesh *mesh = nullptr,
            bool verbose = false,
            bool showStats = false);
  void _centerByDrift(const VectorDouble& dataVect,int ivar=0,bool useSel=true) const;
  void _computeDriftCoeffs() const;
  void _purge();
  bool _isSimulationRequested() const;
  bool _isKrigingRequested() const;
  void _computeLk() const;
  void _computeKriging() const;
  void _computeSimuNonCond() const;
  void _computeSimuCond() const;
  void _addNuggetOnResult(VectorDouble &result) const;
  void _addDrift(Db* db, VectorDouble &result, int ivar = 0, bool useSel = true);
  void _setUseCholesky(int useCholesky = -1, bool verbose = false);
  double _computeLogLikelihood(int nbsimu = 1) const;
  #ifndef SWIG
    static void _projecLocal(Db* dbout,
                             const AMesh* meshing,
                             std::vector<double>& working,
                             VectorDouble& result);
  #endif

private:
  const Db*                    _data; // External Pointer
  ESPDECalcMode                _calcul;
  PrecisionOpMultiConditional* _precisionsKrig;
  PrecisionOpMultiConditional* _precisionsSimu;
  std::vector<PrecisionOp*>    _pilePrecisions; // Dimension: number of valid covariances
  std::vector<ProjMatrix*>     _pileProjMatrix; // Dimension: number of valid covariances
  std::vector<const AMesh*>    _meshingSimu;    // Dimension: number of valid covariances
  std::vector<const AMesh*>    _meshingKrig;    // Dimension: number of valid covariances
  mutable VectorDouble         _driftCoeffs;
  Model*                       _model; // External pointer
  mutable std::vector<std::vector<double>>   _workingKrig;     // Number of Mesh apices * Number of valid covariances
  mutable std::vector<std::vector<double>>   _workingSimu;     // Number of Mesh apices * Number of valid covariances
  mutable std::vector<double>                _workingData;     // Number of valid data
  mutable std::vector<double>                _workingDataInit; // Number of valid data
  std::vector<ProjMatrix*>     _projOnDbOut;
  VectorInt                    _adressesICov;
  double _nugget;
  VectorVectorDouble _driftTab;
  bool _requireCoeffs;
  mutable bool _isCoeffsComputed;
  bool _deleteMesh;
  bool _useCholesky;

  SPDEParam _params;
};

GSTLEARN_EXPORT int krigingSPDE(Db *dbin,
                                Db *dbout,
                                Model *model,
                                bool flag_est = true,
                                bool flag_std = false,
                                const AMesh* mesh = nullptr,
                                int useCholesky = -1,
                                const SPDEParam& params = SPDEParam(),
                                int nbMC = 10,
                                bool verbose = false,
                                bool showStats = false,
                                const NamingConvention &namconv = NamingConvention("KrigingSPDE"));
GSTLEARN_EXPORT int simulateSPDE(Db *dbin,
                                 Db *dbout,
                                 Model *model,
                                 int nbsimu = 1,
                                 const AMesh *mesh = nullptr,
                                 int useCholesky = -1,
                                 const SPDEParam& params = SPDEParam(),
                                 bool verbose = false,
                                 bool showStats = false,
                                 const NamingConvention &namconv = NamingConvention("SimuSPDE"));
GSTLEARN_EXPORT double logLikelihoodSPDE(Db *dbin,
                                         Db *dbout,
                                         Model *model,
                                         const AMesh *mesh = nullptr,
                                         int useCholesky = -1,
                                         int nbsimu = 1,
                                         const SPDEParam& params = SPDEParam(),
                                         bool verbose = false);
GSTLEARN_EXPORT MatrixSparse* buildInvNugget(Db *dbin, Model *model, const SPDEParam& params = SPDEParam());


GSTLEARN_EXPORT VectorDouble krigingSPDENew(
  Db* dbin,
  Db* dbout,
  Model* model,
  bool flag_est = true,
  bool flag_std = false,
  const VectorMeshes& meshes      = VectorMeshes(),
  int useCholesky                 = -1,
  const SPDEParam& params = SPDEParam(),
  int nbMC = 10,
  bool verbose                    = false,
  const NamingConvention& namconv = NamingConvention("KrigingSPDE"));

                             