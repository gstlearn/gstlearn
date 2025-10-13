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
#include "API/SPDE.hpp"
#include "API/SPDEParam.hpp"
#include "Basic/Law.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Enum/ECov.hpp"
#include "LinearOp/InvNuggetOp.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"
#include "LinearOp/PrecisionOpMultiMatrix.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "LinearOp/ProjMultiMatrix.hpp"
#include "LinearOp/SPDEOp.hpp"
#include "LinearOp/SPDEOpMatrix.hpp"
#include "LinearOp/ShiftOpMatrix.hpp"
#include "LithoRule/RuleProp.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Model/Model.hpp"
#include "geoslib_define.h"

#include <cmath>

namespace gstlrn
{

/**
 * The class constructor with the following arguments:
 *
 * @param dbin Input file (nullptr for non conditional simulation)
 * @param model  This compulsory argument is a LMC of Matern's (or Markov?) basic structures with possibly a nugget effect
 * @param flagKrig TRUE for Kriging (or conditional simulation)
 * @param flagSimu TRUE for simulations (either conditional or non-conditional)
 * @param useCholesky Define the choice regarding Cholesky
 * @param params Set of SPDE parameters
 * @param domain Db used to define the domain for meshing (in union with dbin)
 */
SPDE::SPDE(const Db* dbin,
           Model* model,
           bool flagKrig,
           bool flagSimu,
           Id useCholesky,
           const SPDEParam& params,
           const Db* domain)
  : _dbin(dbin)
  , _dbout()
  , _domain(domain)
  , _model(model)
  , _meshesKInit(nullptr)
  , _meshesSInit(nullptr)
  , _projInKInit(nullptr)
  , _projInSInit(nullptr)
  , _projOutKInit(nullptr)
  , _projOutSInit(nullptr)
  , _invnoiseobjInit(nullptr)
  , _isReady(false)
  , _flagCholesky(false)
  , _flagKrig(flagKrig)
  , _flagSimu(flagSimu)
  , _driftCoeffs()
  , _createMeshesK(false)
  , _meshesK()
  , _createMeshesS(false)
  , _meshesS()
  , _createAinK(false)
  , _AinK(nullptr)
  , _createAinS(false)
  , _AinS(nullptr)
  , _createAoutK(false)
  , _AoutK(nullptr)
  , _createAoutS(false)
  , _AoutS(nullptr)
  , _QopK(nullptr)
  , _QopS(nullptr)
  , _Qom(nullptr)
  , _createInvNoise(false)
  , _invnoiseobj(nullptr)
  , _spdeop(nullptr)
  , _params(params)
{
  _defineFlagCholesky(useCholesky, model);
}

SPDE::~SPDE()
{
  _cleanMeshes(true);
  _cleanMeshes(false);

  _cleanProjection(true, false);
  _cleanProjection(true, true);
  _cleanProjection(false, false);
  _cleanProjection(false, true);

  _cleanSpdeOperator();
}

void SPDE::_cleanSpdeOperator()
{
  delete _Qom;
  _Qom = nullptr;
  delete _QopK;
  _QopK = nullptr;
  delete _QopS;
  _QopS = nullptr;
  if (_createInvNoise) delete _invnoiseobj;
  _invnoiseobj = nullptr;
  delete _spdeop;
  _spdeop = nullptr;
}

/**
 * Define if Cholesky must be used or not
 * @param useCholesky: 1 for YES; 0 for No; -1: set optimal default (according to NDim)
 * @param model: The ModelGeneric to use
 * @param verbose: Verbose flag
 */
void SPDE::_defineFlagCholesky(Id useCholesky, const Model* model, bool verbose)
{
  if (useCholesky == -1)
  {
    _flagCholesky = model->getNDim() == 2;
  }
  else
    _flagCholesky = (useCholesky == 1);

  if (verbose)
  {
    if (_flagCholesky)
      message("- Choice for the Cholesky option = ON");
    else
      message("- Choice for the Cholesky option = OFF");
    if (useCholesky == -1)
      message(" (Automatic setting)\n");
    else
      message("\n");
  }
}

Id SPDE::setMeshes(bool flagForKrig, const VectorMeshes* meshes)
{
  if (meshes == nullptr) return 0;
  if (_model == nullptr)
  {
    messerr("The 'model' shoudl be defined beforehand");
    return 1;
  }
  Id ncov  = _model->getNCov(true);
  Id nmesh = static_cast<Id>(meshes->size());

  // Preliminary check
  if (nmesh != 1 && nmesh != ncov)
  {
    messerr("The number of meshes (%d) should be either 1", nmesh);
    messerr("or equal to the number of covariances (%d) contained in the Model", ncov);
    messerr("(Nugget Effect excluded)");
    return 1;
  }
  if (flagForKrig)
    _meshesKInit = meshes;
  else
    _meshesSInit = meshes;

  _isReady = false;
  return 0;
}

Id SPDE::setProjIn(bool flagForKrig, const ProjMultiMatrix* proj)
{
  if (proj == nullptr) return 0;
  if (flagForKrig)
  {
    if (_meshesKInit == nullptr)
    {
      messerr("You may not define specific Input Projection for Kriging");
      messerr("as you have not defined specific Meshes for Kriging beforehand");
      return 1;
    }
    if (!_isValidProjection(_dbin, _meshesKInit, proj)) return 1;
    _projInKInit = proj;
  }
  else
  {
    if (_meshesSInit == nullptr)
    {
      messerr("You may not define specific Input Projection for Simulations");
      messerr("as you have not defined specific Meshes for Simulations beforehand");
      return 1;
    }
    if (!_isValidProjection(_dbin, _meshesSInit, proj)) return 1;
    _projInSInit = proj;
  }
  _isReady = false;
  return 0;
}

Id SPDE::setInvNoise(const ASimulable* invnoise)
{
  if (invnoise == nullptr) return 1;
  if (_dbin == nullptr) return 1;
  _invnoiseobjInit = invnoise;
  _isReady         = false;
  return 0;
}

/**
 * @brief Define the output Db as well as specific Projection matrices (optional)
 *
 * @param dbout  Pointer of the output 'Db'
 * @param projK  Projection matrix for Kriging (optional)
 * @param projS  Projection matrix for Simulations (optional)
 * @param flagApply see Notes
 * @param verbose  Verbose flag
 * @return Id
 *
 * @note If 'dbout' is not defined, this function does nothing
 * @note 'flagApply' is set to TRUE when this method is used on a 'Ready' SPDE system
 * for defining a new 'dbout' for example.
 */
Id SPDE::setDbAndProjOut(Db* dbout,
                         const ProjMultiMatrix* projK,
                         const ProjMultiMatrix* projS,
                         bool flagApply,
                         bool verbose)
{
  if (dbout == nullptr) return 0;
  _dbout = dbout;

  // For Kriging
  if (projK != nullptr)
  {
    if (_meshesKInit == nullptr)
    {
      messerr("You may not define specific Output Projection for Kriging");
      messerr("as you have not defined specific Meshes for Kriging beforehand");
      return 1;
    }
    if (!_isValidProjection(_dbout, _meshesKInit, projK)) return 1;
    _projOutKInit = projK;
  }

  // For simulations
  if (projS != nullptr)
  {
    if (_meshesSInit == nullptr)
    {
      messerr("You may not define specific Projection for Simulations");
      messerr("as you have not defined specific Meshes for Simulations beforehand");
      return 1;
    }
    if (!_isValidProjection(_dbout, _meshesSInit, projS)) return 1;
    _projOutSInit = projS;
  }

  // Apply the output projection as in makeReady function
  if (flagApply)
  {
    if (_flagKrig)
    {
      if (_defineProjection(false, true, verbose)) return 1;
    }

    if (_flagSimu)
    {
      if (_defineProjection(false, false, verbose)) return 1;
    }
  }

  return 0;
}

bool SPDE::_isValidProjection(const Db* db, const VectorMeshes* meshes, const ProjMultiMatrix* proj)
{
  // Get the number of structures and variables
  Id ncov = _model->getNCov(true);
  Id nvar = _model->getNVar();

  // The projection matrix is provided: check consistency of its dimensions
  Id npoint  = proj->getNPoint();
  Id napices = proj->getNApex();

  Id ncol = 0;
  for (Id icov = 0; icov < ncov; icov++)
  {
    ncol += nvar * (*meshes)[icov]->getNApices();
  }
  if (ncol != napices)
  {
    messerr("Number of Columns of 'proj' (%d) is not correct", napices);
    messerr("- Number of covariances = %d", ncov);
    messerr("- Number of variables = %d", nvar);
    for (Id icov = 0; icov < nvar; icov++)
      messerr("- Number of apices for Meshing (%ld) = %d",
              icov + 1, (*meshes)[icov]->getNApices());
    return false;
  }

  Id nrow = 0;
  for (Id ivar = 0; ivar < nvar; ivar++)
    nrow += db->getNSampleActiveAndDefined(ivar);
  if (nrow != npoint)
  {
    messerr("Number of Rows of 'projm' (%d) is not correct", npoint);
    messerr("- Number of variables = %d", nvar);
    for (Id ivar = 0; ivar < nvar; ivar++)
      messerr("- Number of samples for variable %ld = %l\n",
              ivar, db->getNSampleActiveAndDefined(ivar));
    return false;
  }
  return true;
}

VectorMeshes SPDE::_duplicateMeshes(bool flagForKrige)
{
  const auto* src = flagForKrige ? _meshesKInit : _meshesSInit;

  Id ncov = _model->getNCov(true);
  VectorMeshes meshes(ncov);
  for (Id icov = 0; icov < ncov; icov++)
    meshes[icov] = (*src)[0];
  return meshes;
}

VectorMeshes SPDE::_defineMeshesFromDbs(bool flagKrige)
{
  // Create the domain (by merging dbin and dbout)
  bool isBuilt          = false;
  const Db* localDomain = Db::coverSeveralDbs(_dbin, _domain, &isBuilt);

  Id refine  = (flagKrige) ? _params.getRefineK() : _params.getRefineS();
  Id ncovtot = _model->getNCov(false);
  Id ncov    = _model->getNCov(true);
  VectorMeshes meshes(ncov);

  // Create as many meshes as they are structures (nugget excluded)
  for (Id jcov = 0, icov = 0; jcov < ncovtot; jcov++)
  {
    if (_model->getCovType(jcov) == ECov::NUGGET) continue;
    const CovAniso* cova = _model->getCovAniso(jcov);
    meshes[icov++]       = MeshETurbo::createFromCova(*cova, localDomain, refine,
                                                      _params.getBorder(),
                                                      _params.isPolarized(), true,
                                                      _params.getNxMax());
  }
  if (isBuilt) delete localDomain;

  return meshes;
}

void SPDE::_printMeshesDetails(const VectorMeshes& meshes)
{
  Id nmesh = static_cast<Id>(meshes.size());
  for (Id imesh = 0; imesh < nmesh; imesh++)
  {
    message("- Mesh #%d/%d: NMeshes = %d, NApices = %d\n",
            imesh + 1, nmesh, meshes[imesh]->getNMeshes(), meshes[imesh]->getNApices());
  }
}

void SPDE::_printProjectionDetails(const ProjMultiMatrix* proj)
{
  if (proj == nullptr) return;
  message("- Number of variables = %d\n", proj->getNVariable());
  message("- Number of Latents   = %d\n", proj->getNLatent());
  message("- Number of Apices    = %d\n", proj->getNApex());
  message("- Number of Points    = %d\n", proj->getNPoint());
}

void SPDE::_cleanMeshes(bool flagForKrige)
{
  if (flagForKrige)
  {
    if (_createMeshesK && !_meshesK.empty())
    {
      for (Id i = 0, n = static_cast<Id>(_meshesK.size()); i < n; i++)
        delete _meshesK[i];
    }
  }
  else
  {
    if (_createMeshesS && !_meshesS.empty())
    {
      for (Id i = 0, n = static_cast<Id>(_meshesS.size()); i < n; i++)
        delete _meshesS[i];
    }
  }
}

void SPDE::_defineMeshes(bool flagForKrige, bool verbose)
{
  const auto* src  = flagForKrige ? _meshesKInit : _meshesSInit;
  auto& dest       = flagForKrige ? _meshesK : _meshesS;
  auto& createMesh = flagForKrige ? _createMeshesK : _createMeshesS;

  Id ncov  = _model->getNCov(true);
  Id nmesh = (src != nullptr) ? static_cast<Id>(src->size()) : 0;

  // Cleaning already existing meshes inforation
  _cleanMeshes(flagForKrige);

  // Particular case of Simulations: if the corresponding mesh is not defined
  // the meshing dedicated to Kriging is used instead.
  if (!flagForKrige && _meshesS.empty() && !_meshesK.empty())
  {
    dest       = _meshesK;
    createMesh = false;
    if (verbose) message("Copy Meshings from Kriging meshings\n");
    return;
  }

  // If the dimensions 'nmesh' and 'ncov' match: copy pointer of the input VectorMeshes
  if (nmesh == ncov)
  {
    dest       = *src;
    createMesh = false;
    if (verbose)
    {
      message("Input Mesh is provided and used as is\n");
      _printMeshesDetails(dest);
    }
    return;
  }

  // Particular case of a single mesh: simply duplicate it (without creating new contents)
  if (nmesh == 1 && nmesh != ncov)
  {
    dest       = _duplicateMeshes(flagForKrige);
    createMesh = false;
    if (verbose)
    {
      message("Duplicating the Input Mesh for each one of the %d covariance(s)\n", ncov);
      _printMeshesDetails(dest);
    }
    return;
  }

  // Generate the Meshes covering the Dbs
  dest       = _defineMeshesFromDbs(flagForKrige);
  createMesh = true;
  if (verbose)
  {
    message("Creating Meshes covering the available Db(s)\n");
    _printMeshesDetails(dest);
  }
}

void SPDE::_cleanProjection(bool flagIn, bool flagForKrige)
{
  const ProjMultiMatrix* dest =
    flagIn
      ? (flagForKrige ? _AinK : _AinS)
      : (flagForKrige ? _AoutK : _AoutS);
  auto& createProj =
    flagIn
      ? (flagForKrige ? _createAinK : _createAinS)
      : (flagForKrige ? _createAoutK : _createAoutS);

  if (createProj) delete dest;
  dest       = nullptr;
  createProj = false;
}

Id SPDE::_defineProjection(bool flagIn, bool flagForKrige, bool verbose)
{
  const auto* db = flagIn ? _dbin : _dbout;
  if (db == nullptr) return 0;

  auto& meshes = flagForKrige ? _meshesK : _meshesS;
  const ProjMultiMatrix* src =
    flagIn
      ? (flagForKrige ? _projInKInit : _projInSInit)
      : (flagForKrige ? _projOutKInit : _projOutSInit);
  const ProjMultiMatrix** dest =
    flagIn
      ? (flagForKrige ? &_AinK : &_AinS)
      : (flagForKrige ? &_AoutK : &_AoutS);
  auto& createProj =
    flagIn
      ? (flagForKrige ? _createAinK : _createAinS)
      : (flagForKrige ? _createAoutK : _createAoutS);

  const String option = flagForKrige ? "Kriging" : "Simulations";
  const String file   = flagIn ? "'dbin'" : "'dbout'";

  Id ncov = _model->getNCov(true);
  Id nvar = _model->getNVar();
  if (verbose) message("Projection of %s for %s: ", file.c_str(), option.c_str());

  // Clean existing projection (if necessary)
  _cleanProjection(flagIn, flagForKrige);

  // Create the new proejction
  if (src == nullptr)
  {
    *dest      = ProjMultiMatrix::createFromDbAndMeshes(db, meshes, ncov, nvar, flagIn);
    createProj = true;
    if (verbose)
    {
      message("Created from Db and Mesh(es)\n");
      _printProjectionDetails(*dest);
    }
  }
  else
  {
    *dest = src;
    if (verbose) message("Copied the input argument\n");
  }

  return 0;
}

Id SPDE::defineSpdeOperator(bool verbose)
{
  // Clean previous material (if necessary)
  _cleanSpdeOperator();

  if (_invnoiseobjInit != nullptr)
  {
    _invnoiseobj    = _invnoiseobjInit;
    _createInvNoise = false;
  }
  else
  {
    _invnoiseobj    = new InvNuggetOp(_dbin, _model, _params, !_flagCholesky);
    _createInvNoise = true;
  }

  if (_flagCholesky)
  {
    VectorMeshes& meshes = _meshesK;
    if (meshes.empty()) meshes = _meshesS;
    _Qom                 = new PrecisionOpMultiMatrix(_model, meshes);
    const auto* invnoise = dynamic_cast<const InvNuggetOp*>(_invnoiseobj);
    if (invnoise == nullptr)
    {
      messerr("You must provide 'invnoise' as a InvNuggetOp for Cholesky case");
      return 1;
    }
    _spdeop = new SPDEOpMatrix(_Qom, _AinK, invnoise);
  }
  else
  {
    _QopK = new PrecisionOpMulti(_model, _meshesK, _params.getUseStencil());
    if (_flagSimu)
      _QopS = new PrecisionOpMulti(_model, _meshesS, _params.getUseStencil());

    _spdeop = new SPDEOp(_QopK, _AinK, _invnoiseobj, _QopS, _AinS);
    _spdeop->setMaxIterations(_params.getCGparams().getNIterMax());
    _spdeop->setTolerance(_params.getCGparams().getEps());
  }

  _spdeop->setVerbose(verbose);
  return 0;
}

Id SPDE::centerDataByDriftInPlace(VectorDouble& Z, bool verbose)
{
  if (Z.empty()) return 1;
  _driftCoeffs.clear();
  bool mustEvaluateDrift = _model->getNDrift() > 0;

  if (mustEvaluateDrift)
  {
    MatrixDense driftMat = _model->evalDriftMatByRanks(_dbin);
    _driftCoeffs         = _spdeop->computeDriftCoeffs(Z, driftMat);
    ASPDEOp::centerDataByDriftMat(Z, driftMat, _driftCoeffs);
    if (verbose) VH::dump("Drift coefficients = ", _driftCoeffs);
  }
  else
  {
    VectorDouble meanVec = _model->evalMeanVecByRanks(_dbin);
    ASPDEOp::centerDataByMeanVec(Z, meanVec);
  }
  return 0;
}

void SPDE::uncenterResultByDriftInPlace(VectorDouble& result) const
{
  Id nvar                = _model->getNVar();
  Id nech                = _dbout->getNSample();
  Id nechred             = _dbout->getNSample(true);
  bool mustEvaluateDrift = (_model->getNDrift() > 0);

  // Loop on the samples of the output file

  double value;
  Id ncols = 0;
  MatrixDense mat;
  for (Id jech = 0, iech = 0; jech < nech; jech++)
  {
    if (!_dbout->isActive(jech)) continue;

    if (mustEvaluateDrift)
    {
      if (_driftCoeffs.empty()) continue;
      (void)_model->evalDriftMatByTargetInPlace(mat, _dbout, jech);
      ncols = mat.getNCols();
    }

    // Loop on the variables

    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      value = 0.;
      if (mustEvaluateDrift)
      {
        for (Id icol = 0; icol < ncols; icol++)
          value += mat.getValue(ivar, icol) * _driftCoeffs[icol];
      }
      else
        value = _model->getMean(ivar);
      result[ivar * nechred + iech] += value;
    }
    iech++;
  }
}

/**
 * @brief Simulate nugget component and add it to one multivariate simulation
 *
 * @param result The array containing ONE multivariate simulation
 *
 * @note The nugget component is added to each variable with its correct variance
 * @note but not accounting for the dependency across different variables.
 */
void SPDE::addNuggetToResult(VectorDouble& result) const
{
  Id rankNugget = _model->getRankNugget();
  if (rankNugget < 0) return;
  Id nvar = _model->getNVar();
  for (Id ivar = 0, ecr = 0; ivar < nvar; ivar++)
  {
    double nugget = _model->getSill(rankNugget, ivar, ivar);
    for (Id iech = 0, nech = static_cast<Id>(result.size()); iech < nech; iech++, ecr++)
      result[iech] += law_gaussian(0., sqrt(nugget));
  }
}

Id SPDE::makeReady(bool verbose)
{
  // Define Meshes
  if (_flagKrig)
  {
    _defineMeshes(true, verbose);
  }
  if (_flagSimu)
  {
    _defineMeshes(false, verbose);
  }

  // Define the projections on Dbin
  if (_flagKrig)
  {
    if (_defineProjection(true, true, verbose)) return 1;
  }
  if (_flagSimu)
  {
    if (_defineProjection(true, false, verbose)) return 1;
  }

  // Define the projections on Dbout
  if (_flagKrig)
  {
    if (_defineProjection(false, true, verbose)) return 1;
  }
  if (_flagSimu)
  {
    if (_defineProjection(false, false, verbose)) return 1;
  }

  // Define the spde operator
  if (defineSpdeOperator(verbose)) return 1;

  _isReady = true;
  return 0;
}

/**
 * Derive the global trend in the SPDE framework
 *
 * @param dbin Input Db (must contain the variable to be estimated)
 * @param model Model definition
 * @param useCholesky Define the choice regarding Cholesky (see _defineCholesky)
 * @param meshesK Meshes description (optional)
 * @param projInK Matrix of projection (optional)
 * @param params Set of SPDE parameters
 * @param verbose Verbose flag
 *
 * @return Returned vector
 */
VectorDouble trendSPDE(Db* dbin,
                       Model* model,
                       Id useCholesky,
                       const VectorMeshes* meshesK,
                       const ProjMultiMatrix* projInK,
                       const SPDEParam& params,
                       bool verbose)
{
  // Preliminary checks
  if (dbin == nullptr) return VectorDouble();
  if (model == nullptr) return VectorDouble();

  // Instantiate SPDE class
  SPDE spde(dbin, model, true, false, useCholesky, params);
  spde.setMeshes(true, meshesK);
  spde.setProjIn(true, projInK);
  if (verbose) mestitle(1, "Trend in SPDE framework (Cholesky=%ld)", static_cast<Id>(spde.getFlagCholesky()));

  // Define Meshes
  if (spde.makeReady(verbose)) return VectorDouble();

  // Read information from the input Db and center it
  VectorDouble Z = spde.getData(true, verbose);
  if (Z.empty()) return VectorDouble();

  return spde.getDriftCoefficients();
}

/**
 * Perform the estimation by KRIGING under the SPDE framework
 *
 * @param dbin Input Db (must contain the variable to be estimated)
 * @param dbout Output Db where the estimation must be performed
 * @param model Model definition
 * @param flag_est True for the estimation
 * @param flag_std True for the standard deviation of estimation error
 * @param useCholesky Define the choice regarding Cholesky (see _defineCholesky)
 * @param meshesK Meshes description (optional)
 * @param projInK Matrix of projection (optional)
 * @param meshesS Meshes used for Variance calulcation (optional)
 * @param projInS Matrix of projection used for Variance calculation (optional)
 * @param projOutK Matrix of projection on dbout used for Kriging (optional)
 * @param projOutS Matrix of projection on dbout used for Simulations (optional)
 * @param params Set of SPDE parameters
 * @param verbose Verbose flag
 * @param namconv Naming convention
 *
 * @return Returned vector
 *
 * @remark Algorithm for 'meshesK' or 'mesheS' and 'projInK' or 'projInS':
 * - Each one of the previous arguments is considered individually and sequentially.
 * - For 'meshes' in general ('meshesK' or 'meshesS'):
 *   - If it is not defined, it is created from 'model' and so as to cover both 'dbin' and 'dbout' (if provided).
 *   - If is already exist, the number of meshes must be equal:
 *     - either to 1: the same mesh is used for all structures
 *     - or to the number of structures (nugget discarded)
 *   - Otherwise an error is raised
 * - For 'projIn' in general ('projInK' or 'projInS'):
 *   - If it is not defined, it is created from 'meshesK'
 *   - If it is already exist, the number of projectors must be equal to the number of meshes
 *   - Otherwise an error is raised
 * - If 'mesheS' does not exist, 'meshesK' is used instead
 * - If 'projInS' does not exist, 'projInK' is used instead
 *
 * 'meshesS' and 'projInS' are used only if 'flag_std' is true and useCholesky=0
 */
Id krigingSPDE(Db* dbin,
               Db* dbout,
               Model* model,
               bool flag_est,
               bool flag_std,
               Id useCholesky,
               const VectorMeshes* meshesK,
               const ProjMultiMatrix* projInK,
               const VectorMeshes* meshesS,
               const ProjMultiMatrix* projInS,
               const ProjMultiMatrix* projOutK,
               const ProjMultiMatrix* projOutS,
               const SPDEParam& params,
               bool verbose,
               const NamingConvention& namconv)
{
  // Preliminary checks
  if (dbin == nullptr) return 1;
  if (dbout == nullptr) return 1;
  if (model == nullptr) return 1;

  // Instantiate SPDE class
  SPDE spde(dbin, model, true, flag_std, useCholesky, params, dbout);
  spde.setMeshes(true, meshesK);
  spde.setMeshes(false, meshesS);
  spde.setProjIn(true, projInK);
  spde.setProjIn(false, projInS);
  spde.setDbAndProjOut(dbout, projOutK, projOutS, false);
  if (verbose) mestitle(1, "Kriging in SPDE framework (Cholesky=%d)", static_cast<Id>(spde.getFlagCholesky()));

  // Set the SPDE class ready
  if (spde.makeReady(verbose)) return 1;

  // Local variables
  Id nvar = model->getNVar();

  // Read information from the input Db and center it
  VectorDouble Z = spde.getData(true, verbose);
  if (Z.empty()) return 1;

  // Performing the task and storing results in 'dbout'
  // This is performed in ONE step to avoid additional core allocation

  if (flag_est)
  {
    VectorDouble result = spde.kriging(Z);
    Id iuid             = dbout->addColumns(result, "estim", ELoc::Z, 0, true, 0., nvar);
    namconv.setNamesAndLocators(dbin, VectorString(), ELoc::Z, nvar, dbout, iuid,
                                "estim");
  }
  if (flag_std)
  {
    VectorDouble result = spde.stdev(Z);
    Id iuid             = dbout->addColumns(result, "stdev", ELoc::UNKNOWN, 0, true, 0., nvar);
    namconv.setNamesAndLocators(dbin, VectorString(), ELoc::Z, nvar, dbout, iuid,
                                "stdev");
  }

  return 0;
}

/**
 * Perform the conditional SIMULATIONs in the SPDE framework
 *
 * @param dbin Input Db (variable to be estimated). Only for conditional simulations
 * @param dbout Output Db where the estimation must be performed
 * @param model Model definition
 * @param nbsimu Number of simulations
 * @param useCholesky Define the choice regarding Cholesky (see _defineCholesky)
 * @param meshesK Meshes used for Kriging (optional)
 * @param projInK Matrix of projection used for Kriging (optional)
 * @param meshesS Meshes used for Simulations (optional)
 * @param projInS Matrix of projection used for Simulations (optional)
 * @param projOutK Matrix of projection on dbout used for Kriging (optional)
 * @param projOutS Matrix of projection on dbout used for Simulations (optional)
 * @param params Set of SPDE parameters
 * @param verbose Verbose flag
 * @param namconv see NamingConvention
 *
 * @return Error returned code
 *
 * @remark Algorithm for 'meshesK', 'projInK', 'meshesS' and 'projInS':
 * - Each one of the previous arguments is considered individually and sequentially.
 * - For 'meshes' in general ('meshesK' or 'meshesS'):
 *   - If it is not defined, it is created from 'model' and so as to cover both 'dbin' and 'dbout' (if provided).
 *   - If is already exist, the number of meshes must be equal:
 *     - either to 1: the same mesh is used for all structures
 *     - or to the number of structures (nugget discarded)
 *   - Otherwise an error is raised
 * - For 'projIn' in general ('projInK' or 'projInS'):
 *   - If it is not defined, it is created from 'meshes'
 *   - If it is already exist, the number of projectors must be equal to the number of meshes
 *   - Otherwise an error is raised
 * - If 'mesheS' does not exist, 'meshesK' is used instead
 * - If 'projInS' does not exist, 'projInK' is used instead
 */
Id simulateSPDE(Db* dbin,
                Db* dbout,
                Model* model,
                Id nbsimu,
                Id useCholesky,
                const VectorMeshes* meshesK,
                const ProjMultiMatrix* projInK,
                const VectorMeshes* meshesS,
                const ProjMultiMatrix* projInS,
                const ProjMultiMatrix* projOutK,
                const ProjMultiMatrix* projOutS,
                const SPDEParam& params,
                bool verbose,
                const NamingConvention& namconv)
{
  if (dbout == nullptr) return 1;
  if (model == nullptr) return 1;

  // Instantiate SPDE class
  bool flagKrig = dbin != nullptr;
  SPDE spde(dbin, model, flagKrig, true, useCholesky, params, dbout);
  spde.setMeshes(true, meshesK);
  spde.setMeshes(false, meshesS);
  spde.setProjIn(true, projInK);
  spde.setProjIn(false, projInS);
  spde.setDbAndProjOut(dbout, projOutK, projOutS, false);
  if (verbose) mestitle(1, "Simulation in SPDE framework (Conditional=%d, Cholesky=%d)",
                        dbin != nullptr, static_cast<Id>(spde.getFlagCholesky()));

  // Set the SPDE class ready
  if (spde.makeReady(verbose)) return 1;

  // Local variables
  Id nvar    = model->getNVar();
  Id nechred = dbout->getNSample(true);

  VectorDouble Z;
  VectorDouble driftCoeffs;
  if (spde.getFlagKrig())
  {
    // Read information from the input Db and center it
    Z = spde.getData(true, verbose);
    if (Z.empty()) return 1;
  }

  // Perform the Simulation and storage.
  // All is done in ONE step to avoid additional storage

  Id iuid = dbout->addColumnsByConstant(nvar * nbsimu);
  VectorDouble local(nechred);
  VectorDouble result;

  for (Id isimu = 0; isimu < nbsimu; isimu++)
  {
    result = spde.simulate(Z);

    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      VH::extractInPlace(result, local, ivar * nechred);
      Id juid = iuid + ivar * nbsimu + isimu;
      dbout->setColumnByUID(local, juid, true);
    }
  }
  namconv.setNamesAndLocators(dbin, VectorString(), ELoc::Z, nvar, dbout, iuid,
                              "", nbsimu);

  return 0;
}

/**
 * Perform the conditional SIMULATIONs using PGS technique in the SPDE framework
 *
 * @param dbin Input Db (variable to be estimated). Only for conditional simulations
 * @param dbout Output Db where the estimation must be performed
 * @param model Model definition
 * @param ruleprop RuleProp structure describing the Rule and the Proportions
 * @param nbsimu Number of simulations
 * @param useCholesky Define the choice regarding Cholesky (see _defineCholesky)
 * @param meshesK Meshes used for Kriging (optional)
 * @param projInK Matrix of projection used for Kriging (optional)
 * @param meshesS Meshes used for Simulations (optional)
 * @param projInS Matrix of projection used for Simulations (optional)
 * @param projOutK Matrix of projection on dbout used for Kriging (optional)
 * @param projOutS Matrix of projection on dbout used for Simulations (optional)
 * @param params Set of SPDE parameters
 * @param verbose Verbose flag
 * @param namconv see NamingConvention
 *
 * @return Error returned code
 *
 * @remark Algorithm for 'meshesK', 'projInK', 'meshesS' and 'projInS':
 * - Each one of the previous arguments is considered individually and sequentially.
 * - For 'meshes' in general ('meshesK' or 'meshesS'):
 *   - If it is not defined, it is created from 'model' and so as to cover both 'dbin' and 'dbout' (if provided).
 *   - If is already exist, the number of meshes must be equal:
 *     - either to 1: the same mesh is used for all structures
 *     - or to the number of structures (nugget discarded)
 *   - Otherwise an error is raised
 * - For 'projIn' in general ('projInK' or 'projInS'):
 *   - If it is not defined, it is created from 'meshes'
 *   - If it is already exist, the number of projectors must be equal to the number of meshes
 *   - Otherwise an error is raised
 * - If 'mesheS' does not exist, 'meshesK' is used instead
 * - If 'projInS' does not exist, 'projInK' is used instead
 *
 * @remark Note that, at this stage of development of this method, the data must be provided
 * in Gaussian scale.
 */
Id simPGSSPDE(Db* dbin,
              Db* dbout,
              Model* model,
              const RuleProp& ruleprop,
              Id nbsimu,
              Id useCholesky,
              const VectorMeshes* meshesK,
              const ProjMultiMatrix* projInK,
              const VectorMeshes* meshesS,
              const ProjMultiMatrix* projInS,
              const ProjMultiMatrix* projOutK,
              const ProjMultiMatrix* projOutS,
              const SPDEParam& params,
              bool verbose,
              const NamingConvention& namconv)
{
  if (dbout == nullptr) return 1;
  if (model == nullptr) return 1;

  bool flagKrig = dbin != nullptr;
  SPDE spde(dbin, model, flagKrig, true, useCholesky, params, dbout);
  spde.setMeshes(true, meshesK);
  spde.setMeshes(false, meshesS);
  spde.setProjIn(true, projInK);
  spde.setProjIn(false, projInS);
  spde.setDbAndProjOut(dbout, projOutK, projOutS, false);
  if (verbose) mestitle(1, "PluriGaussian Simulation in SPDE framework (Conditional=%d, Cholesky=%d)",
                        dbin != nullptr, static_cast<Id>(spde.getFlagCholesky()));

  // Set the SPDE class ready
  if (spde.makeReady(verbose)) return 1;

  // Local variables
  Id nvar    = model->getNVar();
  Id nechred = dbout->getNSample(true);

  VectorDouble Z;
  if (spde.getFlagKrig())
  {
    // Read information from the input Db and center it
    Z = spde.getData(true, verbose);
    if (Z.empty()) return 1;
  }

  // Perform the Simulation and storage.
  // All is done in ONE step to avoid additional storage
  Id iuid = dbout->addColumnsByConstant(nvar * nbsimu);
  VectorDouble local(nechred);
  VectorDouble result;

  // Loop on the simulations
  for (Id isimu = 0; isimu < nbsimu; isimu++)
  {
    result = spde.simulate(Z);

    // Loop on the variables
    VectorInt iuids;
    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      VH::extractInPlace(result, local, ivar * nechred);
      Id juid = iuid + ivar * nbsimu + isimu;
      dbout->setColumnByUID(local, juid, true);
      iuids.push_back(juid);
    }

    // Convert the resulting simulation into categories
    dbout->setLocatorsByUID(iuids, ELoc::SIMU);
    ruleprop.gaussToCategory(dbout, namconv);
    dbout->deleteColumnsByUID(iuids);
  }

  return 0;
}

/**
 * Calculate the Log-Likelihood under the SPDE framework
 *
 * @param dbin Input Db (must contain the variable to be estimated)
 * @param model Model definition
 * @param useCholesky Define the choice regarding Cholesky (see _defineCholesky)
 * @param meshes Meshes description (optional)
 * @param projIn Matrix of projection (optional)
 * @param params Set of SPDE parameters
 * @param verbose True for verbose output
 * @return Returned value
 */
double logLikelihoodSPDE(Db* dbin,
                         Model* model,
                         Id useCholesky,
                         const VectorMeshes* meshes,
                         const ProjMultiMatrix* projIn,
                         const SPDEParam& params,
                         bool verbose)
{
  if (dbin == nullptr) return 1;
  if (model == nullptr) return 1;
  if (dbin->getNLoc(ELoc::Z) != 1)
  {
    messerr("'dbin' must contain ONE variable (Z locator)");
    return 1;
  }

  // Instantiate SPDE class
  SPDE spde(dbin, model, true, false, useCholesky, params);
  spde.setMeshes(true, meshes);
  spde.setProjIn(true, projIn);
  if (verbose) mestitle(1, "Log-likelhood calculation in SPDE framework( Cholesky=%d)",
                        static_cast<Id>(spde.getFlagCholesky()));

  // Set the SPDE class ready
  if (spde.makeReady(verbose)) return 1;

  // Read information from the input Db and center it
  VectorDouble Z = spde.getData(true, verbose);
  if (Z.empty()) return 1;

  // Performing the task
  return spde.loglikelihood(Z, verbose);
}

const MatrixSparse* SPDE::getQ() const
{
  if (!_isReady)
  {
    messerr("You must call SPDE::makeReady() beforehand");
    return nullptr;
  }
  auto* spdeopmatrix = dynamic_cast<SPDEOpMatrix*>(_spdeop);
  if (spdeopmatrix == nullptr) return nullptr;

  const auto* mat = dynamic_cast<const PrecisionOpMultiMatrix*>(spdeopmatrix->getQKriging());
  return mat->getQ();
}

const MatrixSparse* SPDE::getInvNoise() const
{
  if (!_isReady)
  {
    messerr("You must call SPDE::makeReady() beforehand");
    return nullptr;
  }
  const auto* invnoise = dynamic_cast<const InvNuggetOp*>(_invnoiseobj);
  if (invnoise == nullptr) return nullptr;

  return invnoise->getInvNuggetMatrix();
}

const MatrixSparse* SPDE::getProj() const
{
  if (!_isReady)
  {
    messerr("You must call SPDE::makeReady() beforehand");
    return nullptr;
  }
  auto* spdeopmatrix = dynamic_cast<SPDEOpMatrix*>(_spdeop);
  if (spdeopmatrix == nullptr) return nullptr;

  const auto* mat = dynamic_cast<const ProjMultiMatrix*>(spdeopmatrix->getProjKriging());
  return mat->getProj();
}

const PrecisionOpMulti* SPDE::getPrecisionKrig() const
{
  if (!_isReady)
  {
    messerr("You must call SPDE::makeReady() beforehand");
    return nullptr;
  }
  return _flagCholesky ? _Qom : _QopK;
}

VectorDouble SPDE::getData(bool flagCenter, bool verbose)
{
  VectorDouble Z;
  if (!_isReady)
  {
    messerr("You must call SPDE::makeReady() beforehand");
    return Z;
  }
  Z = _dbin->getColumnsActiveAndDefined(ELoc::Z);

  if (flagCenter)
  {
    if (centerDataByDriftInPlace(Z, verbose)) return 1;
  }
  return Z;
}

VectorDouble SPDE::kriging(const VectorDouble& Z) const
{
  if (!_isReady)
  {
    messerr("You must call SPDE::makeReady() beforehand");
    return VectorDouble();
  }
  VectorDouble result = getSPDEOp()->kriging(Z, _AoutK);
  uncenterResultByDriftInPlace(result);
  return result;
}

VectorDouble SPDE::stdev(const VectorDouble& Z) const
{
  if (!_isReady)
  {
    messerr("You must call SPDE::makeReady() beforehand");
    return VectorDouble();
  }
  VectorDouble result = getSPDEOp()->stdev(Z, getNMC(), getSeed(),
                                           getAoutK(), getAoutS());
  return result;
}

VectorDouble SPDE::simulate(const VectorDouble& Z) const
{
  VectorDouble result = (getFlagKrig())
                        ? getSPDEOp()->simCond(Z, getAoutK(), getAoutS())
                        : getSPDEOp()->simNonCond(getAoutS());
  addNuggetToResult(result);
  uncenterResultByDriftInPlace(result);
  return result;
}

double SPDE::loglikelihood(const VectorDouble& Z, bool verbose) const
{
  Id size        = static_cast<Id>(Z.size());
  double logdet  = getSPDEOp()->computeTotalLogDet(getNMC(), getSeed());
  double quad    = getSPDEOp()->computeQuadratic(Z);
  double loglike = TEST;
  if (!FFFF(logdet) && !FFFF(quad))
    loglike = -0.5 * (logdet + quad + size * log(2. * GV_PI));

  if (verbose)
  {
    message("Likelihood calculation:\n");
    message("Nb. active samples = %d\n", size);
    message("Nb. Monte-Carlo    = %d\n", getNMC());
    message("Cholesky           = %d\n", getFlagCholesky());
    message("Log-Determinant    = %lf\n", logdet);
    message("Quadratic term     = %lf\n", quad);
    message("Log-likelihood     = %lf\n", loglike);
  }
  return loglike;
}

} // namespace gstlrn
