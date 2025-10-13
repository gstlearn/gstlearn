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

#include "API/SPDEParam.hpp"
#include "Basic/NamingConvention.hpp"
#include "LinearOp/ASimulable.hpp"
#include "LinearOp/InvNuggetOp.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"
#include "LinearOp/PrecisionOpMultiMatrix.hpp"
#include "LinearOp/ProjMultiMatrix.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
class ShiftOpMatrix;
class PrecisionOp;
class Db;
class DbGrid;
class MeshETurbo;
class Model;
class RuleProp;
class SPDEOp;

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
  SPDE(const Db* dbin,
       Model* model,
       bool flagKrig,
       bool flagSimu,
       Id useCholesky          = -1,
       const SPDEParam& params = SPDEParam(),
       const Db* domain        = nullptr);
  SPDE(const SPDE& r)            = delete;
  SPDE& operator=(const SPDE& r) = delete;
  virtual ~SPDE();

public:
  bool getFlagCholesky() const { return _flagCholesky; }
  bool getFlagKrig() const { return _flagKrig; }
  bool getFlagSimu() const { return _flagSimu; }
  const VectorMeshes& getMeshesK() const { return _meshesK; }
  const VectorMeshes& getMeshesS() const { return _meshesS; }
  const ProjMultiMatrix* getAinK() const { return _AinK; }
  const ProjMultiMatrix* getAinS() const { return _AinS; }
  const ProjMultiMatrix* getAoutK() const { return _AoutK; }
  const ProjMultiMatrix* getAoutS() const { return _AoutS; }
  VectorDouble getDriftCoefficients() const { return _driftCoeffs; }

  Id getSeed() const { return _params.getSeedMC(); }
  Id getNMC() const { return _params.getNMC(); }
  const SPDEOp* getSPDEOp() const { return _spdeop; }

  const MatrixSparse* getQ() const;
  const MatrixSparse* getProj() const;
  const PrecisionOpMulti* getPrecisionKrig() const;
  const MatrixSparse* getInvNoise() const;

  Id defineSpdeOperator(bool verbose = false);
  Id centerDataByDriftInPlace(VectorDouble& Z, bool verbose = false);
  void uncenterResultByDriftInPlace(VectorDouble& result) const;
  void addNuggetToResult(VectorDouble& result) const;

  Id setMeshes(bool flagForKrig, const VectorMeshes* meshes);
  Id setProjIn(bool flagForKrig, const ProjMultiMatrix* proj);
  Id setInvNoise(const ASimulable* invnoise);
  Id setDbAndProjOut(Db* dbout,
                     const ProjMultiMatrix* projK = nullptr,
                     const ProjMultiMatrix* projS = nullptr,
                     bool flagApply               = true,
                     bool verbose                 = false);
  Id makeReady(bool verbose = false);
  VectorDouble getData(bool flagCenter, bool verbose = false);
  VectorDouble kriging(const VectorDouble& Z) const;
  VectorDouble stdev(const VectorDouble& Z) const;
  VectorDouble simulate(const VectorDouble& Z) const;
  double loglikelihood(const VectorDouble& Z, bool verbose = false) const;

private:
  void _defineFlagCholesky(Id useCholesky, const Model* model, bool verbose = false);
  VectorMeshes _defineMeshesFromDbs(bool flagKrige);
  VectorMeshes _duplicateMeshes(bool flagForKrige);
  void _defineMeshes(bool flagForKrige, bool verbose);
  Id _defineProjection(bool flagIn, bool flagForKrige, bool verbose);
  static void _printMeshesDetails(const VectorMeshes& meshes);
  static void _printProjectionDetails(const ProjMultiMatrix* proj);
  bool _isValidProjection(const Db* db, const VectorMeshes* meshes, const ProjMultiMatrix* proj);
  void _cleanMeshes(bool flagForKrige);
  void _cleanProjection(bool flagIn, bool flagForKrige);
  void _cleanSpdeOperator();

private:
  // External pointers
  const Db* _dbin;
  const Db* _dbout;
  const Db* _domain;
  Model* _model;
  const VectorMeshes* _meshesKInit;
  const VectorMeshes* _meshesSInit;
  const ProjMultiMatrix* _projInKInit;
  const ProjMultiMatrix* _projInSInit;
  const ProjMultiMatrix* _projOutKInit;
  const ProjMultiMatrix* _projOutSInit;
  const ASimulable* _invnoiseobjInit;

  // Local information
  bool _isReady;
  bool _flagCholesky;
  bool _flagKrig;
  bool _flagSimu;
  VectorDouble _driftCoeffs;
  bool _createMeshesK;
  VectorMeshes _meshesK;
  bool _createMeshesS;
  VectorMeshes _meshesS;
  bool _createAinK;
  const ProjMultiMatrix* _AinK;
  bool _createAinS;
  const ProjMultiMatrix* _AinS;
  bool _createAoutK;
  const ProjMultiMatrix* _AoutK;
  bool _createAoutS;
  const ProjMultiMatrix* _AoutS;
  PrecisionOpMulti* _QopK;
  PrecisionOpMulti* _QopS;
  PrecisionOpMultiMatrix* _Qom;
  bool _createInvNoise;
  const ASimulable* _invnoiseobj;
  SPDEOp* _spdeop;
  SPDEParam _params;
};

GSTLEARN_EXPORT VectorDouble trendSPDE(Db* dbin,
                                       Model* model,
                                       Id useCholesky                 = -1,
                                       const VectorMeshes* meshesK    = nullptr,
                                       const ProjMultiMatrix* projInK = nullptr,
                                       const SPDEParam& params        = SPDEParam(),
                                       bool verbose                   = false);
GSTLEARN_EXPORT Id krigingSPDE(Db* dbin,
                               Db* dbout,
                               Model* model,
                               bool flag_est                   = true,
                               bool flag_std                   = false,
                               Id useCholesky                  = -1,
                               const VectorMeshes* meshesK     = nullptr,
                               const ProjMultiMatrix* projInK  = nullptr,
                               const VectorMeshes* meshesS     = nullptr,
                               const ProjMultiMatrix* projInS  = nullptr,
                               const ProjMultiMatrix* projOutK = nullptr,
                               const ProjMultiMatrix* projOutS = nullptr,
                               const SPDEParam& params         = SPDEParam(),
                               bool verbose                    = false,
                               const NamingConvention& namconv = NamingConvention("KrigingSPDE"));
GSTLEARN_EXPORT Id simulateSPDE(Db* dbin,
                                Db* dbout,
                                Model* model,
                                Id nbsimu                       = 1,
                                Id useCholesky                  = -1,
                                const VectorMeshes* meshesK     = nullptr,
                                const ProjMultiMatrix* projInK  = nullptr,
                                const VectorMeshes* meshesS     = nullptr,
                                const ProjMultiMatrix* projInS  = nullptr,
                                const ProjMultiMatrix* projOutK = nullptr,
                                const ProjMultiMatrix* projOutS = nullptr,
                                const SPDEParam& params         = SPDEParam(),
                                bool verbose                    = false,
                                const NamingConvention& namconv = NamingConvention("SimuSPDE"));
GSTLEARN_EXPORT Id simPGSSPDE(Db* dbin,
                              Db* dbout,
                              Model* model,
                              const RuleProp& ruleprop,
                              Id nbsimu                       = 1,
                              Id useCholesky                  = -1,
                              const VectorMeshes* meshesK     = nullptr,
                              const ProjMultiMatrix* projInK  = nullptr,
                              const VectorMeshes* meshesS     = nullptr,
                              const ProjMultiMatrix* projInS  = nullptr,
                              const ProjMultiMatrix* projOutK = nullptr,
                              const ProjMultiMatrix* projOutS = nullptr,
                              const SPDEParam& params         = SPDEParam(),
                              bool verbose                    = false,
                              const NamingConvention& namconv = NamingConvention("SimPGSSPDE"));
GSTLEARN_EXPORT double logLikelihoodSPDE(Db* dbin,
                                         Model* model,
                                         Id useCholesky                = -1,
                                         const VectorMeshes* meshes    = nullptr,
                                         const ProjMultiMatrix* projIn = nullptr,
                                         const SPDEParam& params       = SPDEParam(),
                                         bool verbose                  = false);

} // namespace gstlrn