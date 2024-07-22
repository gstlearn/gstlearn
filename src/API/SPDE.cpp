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
#include "Enum/ECov.hpp"

#include "API/SPDE.hpp"
#include "Model/ANoStat.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Covariances/CovAniso.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/NamingConvention.hpp"
#include "Model/Model.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "LinearOp/PrecisionOpMultiConditionalCs.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Db/Db.hpp"

#include <math.h>

/**
 * The class constructor with the following arguments:
 *
 * @param model  This compulsory argument is a LMC of Matern's (or Markov?) basic structures with possibly a nugget effect
 * @param domain The Db defining the space dimension and spatial limits where the SPDE model is defined
 * @param data   The Db containing the data for conditioning (optional)
 * @param calcul Option from ESPDECalcMode
 * @param meshUser The mesh for the discretization of the domain
 * @param useCholesky Define the choice regarding Cholesky
 * @param params Set of SPDE parameters
 * @param verbose  Verbose flag
 * @param showStats Display statistics for Linear Operations (when deleting the class)
 *
 * @details
 * Either Domain or a Mesh can be provided:
 * - If a mesh is not provided, the domain is used to define optimal meshes for each structures and
 * for at low/high resolution respectively.
 * (High resolution meshes are used for simulations and low resolution meshes are used for kriging).
 * - If a mesh is provided it is used for all structures and resolutions.
 *
 * The domain or the mesh should have the same spatial reference (ndim, manifold) than the model.
 */
SPDE::SPDE(Model* model,
           const Db* domain,
           const Db* data,
           const ESPDECalcMode& calcul,
           const AMesh* meshUser,
           int useCholesky,
           const SPDEParam& params,
           bool verbose,
           bool showStats)
    : _data(data),
      _calcul(calcul),
      _precisionsKrig(nullptr),
      _precisionsSimu(nullptr),
      _pilePrecisions(),
      _pileProjMatrix(),
      _meshingSimu(),
      _meshingKrig(),
      _driftCoeffs(),
      _model(model),
      _workingKrig(),
      _workingSimu(),
      _workingData(),
      _workingDataInit(),
      _projOnDbOut(),
      _adressesICov(),
      _nugget(0.),
      _driftTab(),
      _requireCoeffs(false),
      _isCoeffsComputed(false),
      _deleteMesh(false),
      _useCholesky(true),
      _params(params)
{
  _setUseCholesky(useCholesky, verbose);

  (void) _init(domain, meshUser, verbose, showStats);
}

SPDE::~SPDE()
{
  _purge();
}

void SPDE::_purge()
{
  delete _precisionsKrig;
  _precisionsKrig = nullptr;
  delete _precisionsSimu;
  _precisionsSimu = nullptr;

  for (int i = 0, n = (int) _pilePrecisions.size(); i < n; i++)
    delete _pilePrecisions[i];
  _pilePrecisions.clear();

  for (int i = 0, n = (int) _pileProjMatrix.size(); i < n; i++)
    delete _pileProjMatrix[i];
  _pileProjMatrix.clear();

  for (int i = 0, n = (int) _projOnDbOut.size(); i < n; i++)
    delete _projOnDbOut[i];
  _projOnDbOut.clear();

  if (_deleteMesh)
  {
    for (int i = 0, n = (int) _meshingSimu.size(); i < n; i++)
      delete _meshingSimu[i];
    for (int i = 0, n = (int) _meshingKrig.size(); i < n; i++)
      delete _meshingKrig[i];
  }
  _meshingSimu.clear();
  _meshingKrig.clear();
}

SPDE* SPDE::create(Model *model,
                   const Db *domain,
                   const Db *data,
                   const ESPDECalcMode &calcul,
                   const AMesh* meshUser,
                   int useCholesky,
                   const SPDEParam& params,
                   bool verbose,
                   bool showStats)
{
  return new SPDE(model, domain, data, calcul, meshUser, useCholesky, params, verbose, showStats);
}

/**
 * Define if Cholesky must be used or not
 * @param useCholesky: 1 for YES; 0 for No; -1: set optimal default
 * @param verbose: Verbose flag
 */
void SPDE::_setUseCholesky(int useCholesky, bool verbose)
{
  if (useCholesky == -1)
  {
    _useCholesky = (_model->getDimensionNumber() == 2);
  }
  else if (useCholesky == 1)
    _useCholesky = true;
  else
    _useCholesky = false;

  // Optional printout
  if (verbose)
  {
    mestitle(1, "SPDE parameters");
    message("- Space dimension = %d\n", _model->getDimensionNumber());
    if (! _meshingKrig.empty())
    {
      for (int imesh = 0; imesh < (int) _meshingKrig.size(); imesh++)
        message("- Number of meshes (Kriging #%d) = %d\n",
                imesh+1, _meshingKrig[imesh]->getNMeshes());
    }
    if (! _meshingSimu.empty())
    {
      for (int imesh = 0; imesh < (int) _meshingSimu.size(); imesh++)
        message("- Number of meshes (Kriging #%d) = %d\n",
                imesh+1, _meshingSimu[imesh]->getNMeshes());
    }

    if (_useCholesky)
      message("- Choice for the Cholesky option = ON");
    else
      message("- Choice for the Cholesky option = OFF");
    if (useCholesky == -1)
      message(" (Automatic setting)\n");
    else
      message("\n");
  }
}

int SPDE::_init(const Db *domain, const AMesh *meshUser, bool verbose, bool showStats)
{
  const ANoStat* nostat = nullptr;
  if (_model->isNoStat()) nostat = _model->getNoStat();

  if (_isKrigingRequested() && _data == nullptr)
  {
    messerr("You must define 'data' when performing Kriging or Conditional Simulations");
    return 1;
  }

  bool useSel = true;
  VectorDouble varianceData;
  double totalSill = 0.;
  PrecisionOp* precision;
  const AMesh* mesh = meshUser;
  ProjMatrix* proj;
  if (_data != nullptr) _driftTab = _model->getDrifts(_data, useSel);
  _requireCoeffs = _driftTab.size() > 0 && _data != nullptr;

  // Allocate the global structures (pointers)
  if (_isSimulationRequested())
  {
    _precisionsSimu = new PrecisionOpMultiConditional();
    _precisionsSimu->mustShowStats(showStats);
  }
  if (_isKrigingRequested() || _requireCoeffs)
  {
    if (_useCholesky)
      _precisionsKrig = new PrecisionOpMultiConditionalCs();
    else
      _precisionsKrig = new PrecisionOpMultiConditional();
    _precisionsKrig->mustShowStats(showStats);
  }

  // Loop on the basic structures
  for (int icov = 0, ncov = _model->getCovaNumber(); icov < ncov; icov++)
  {
    const CovAniso* cova = _model->getCova(icov);
    double sill = cova->getSill(0,0);
    bool flagNoStatRot = false;
    if (nostat != nullptr) flagNoStatRot = nostat->isDefinedforAnisotropy(icov);

    if (cova->getType() == ECov::NUGGET)
    {
      _nugget = sill;
    }
    else if (cova->getType() == ECov::BESSEL_K || cova->getType() == ECov::MARKOV)
    {
      totalSill += sill;

      if (_isSimulationRequested())
      {
        if (meshUser == nullptr)
        {
          mesh = MeshETurbo::createFromCova(*cova, domain, _params.getRefineS(), _params.getBorder(),
                                            useSel, flagNoStatRot, verbose);
          _deleteMesh = true;
        }
        _meshingSimu.push_back(mesh);

        if (_useCholesky)
          precision = new PrecisionOpCs(mesh, _model, icov, false, _params.getCGparams(), verbose);
        else
          precision = new PrecisionOp(mesh, _model, icov, _params.getCGparams(), verbose);
        precision->mustShowStats(showStats);
        _pilePrecisions.push_back(precision);

        proj = new ProjMatrix(_data, mesh, 0);
        _pileProjMatrix.push_back(proj);

        if (_precisionsSimu->push_back(precision, proj) != 0) return 1;
        _precisionsSimu->setVarianceDataVector(varianceData);
        _workingSimu.push_back(VectorDouble(precision->getSize()));
      }

      if (_isKrigingRequested() || _requireCoeffs)
      {
        if (meshUser == nullptr)
        {
          mesh = MeshETurbo::createFromCova(*cova, domain, _params.getRefineK(),
                                            _params.getBorder(), useSel, flagNoStatRot, verbose);
          _deleteMesh = true;
        }
        _meshingKrig.push_back(mesh);

        if (_useCholesky)
          precision = new PrecisionOpCs(mesh, _model, icov, false, _params.getCGparams(), verbose);
        else
          precision = new PrecisionOp(mesh, _model, icov, _params.getCGparams(), verbose);
        precision->mustShowStats(showStats);
        _pilePrecisions.push_back(precision);

        proj = new ProjMatrix(_data, mesh, 0);
        _pileProjMatrix.push_back(proj);

        if (_precisionsKrig->push_back(precision, proj) != 0) return 1;
        _workingKrig.push_back(VectorDouble(precision->getSize()));
      }
    }
    else
    {
      messerr("SPDE is only implemented for Matérn (BESSEL_K) and Markov (MARKOV) covariances");
      return 1;
    }
  }

  // Evaluation of the variance at data point
  if (_isKrigingRequested() && _data != nullptr)
  {
    if (_data->getLocNumber(ELoc::V) > 0)
    {
      // If a variance of measurement error is defined
      // we must intersect it with the definition of the Z-value
      VectorDouble valData = _data->getColumnByLocator(ELoc::Z,0,useSel);
      VectorDouble varData = _data->getColumnByLocator(ELoc::V,0,useSel);
      double eps_loc = _params.getEpsNugget() * totalSill;
      for (int iech = 0; iech < _data->getSampleNumber(true); iech++)
      {
        if (FFFF(valData[iech])) continue;
        {
          double loc_value;
          if (FFFF(varData[iech]))
            loc_value = eps_loc;
          else
          {
            loc_value = MAX(varData[iech], eps_loc);
          }
          varianceData.push_back(loc_value);
        }
      }
    }
    else
    {
      VH::fill(varianceData, MAX(_nugget, _params.getEpsNugget() * totalSill),
               _data->getNumberActiveAndDefined(0));
    }
    _precisionsKrig->setVarianceDataVector(varianceData);

    if (_isSimulationRequested())
      _precisionsSimu->setVarianceDataVector(varianceData);
  }

  return 0;
}

void SPDE::_computeLk() const
{
  VectorVectorDouble rhs = _precisionsKrig->computeRhs(_workingData);
  _precisionsKrig->initLk(rhs, _workingKrig); // Same as evalInverse but with just one iteration
}

void SPDE::_computeKriging() const
{
  VectorVectorDouble rhs = _precisionsKrig->computeRhs(_workingData);
  _precisionsKrig->evalInverse(rhs, _workingKrig);
}

/**
 * Perform one non-conditional simulation on internal meshing
 * The results (for each covariance item) are stored in _workingSimu
 */
void SPDE::_computeSimuNonCond() const
{
  _precisionsSimu->simulateOnMeshings(_workingSimu);
}

/**
 * Perform one conditional simulation on internal meshing
 * The results (for each covariance item) are stored in _workingSimu'
 */
void SPDE::_computeSimuCond() const
{
  // Perform the non conditional simulation on target
  _computeSimuNonCond();

  // Perform the non conditional simulation on data
  VectorDouble temp_dat(_data->getSampleNumber(true));
  _precisionsSimu->simulateOnDataPointFromMeshings(_workingSimu, temp_dat);

  // Calculate the simulation error
  _workingData = _workingDataInit;
  VH::multiplyConstant(temp_dat, -1.);
  VH::addInPlace(_workingData, temp_dat);

  // Conditional Kriging
  _computeKriging();
}

void SPDE::_centerByDrift(const VectorDouble& dataVect,int ivar,bool useSel) const
{
  _computeDriftCoeffs();

  if (_driftCoeffs.empty())
  {
    if (_workingDataInit.empty())
    {
      _workingDataInit.resize(dataVect.size());
    }

    for(int iech = 0, nech = (int) _workingDataInit.size(); iech<nech;iech++)
    {
      _workingDataInit[iech] = dataVect[iech];
    }
  }
  else
  {
    _workingDataInit = _model->evalDriftVarCoefs(_data,_driftCoeffs,ivar,useSel);

    for(int iech = 0, nech = (int) _workingDataInit.size(); iech<nech; iech++)
    {
      _workingDataInit[iech] = dataVect[iech] - _workingDataInit[iech];
    }
  }
}

void SPDE::_addDrift(Db* db, VectorDouble &result, int ivar, bool useSel)
{
  if (! _requireCoeffs) return;
  VectorDouble temp_out = _model->evalDriftVarCoefs(db, _driftCoeffs, ivar, useSel);
  VH::addInPlace(result, temp_out);
}

int SPDE::compute(Db *dbout,
                  int nbsimu,
                  int seed,
                  const NamingConvention &namconv)
{
  VectorDouble dataVect;
  bool useSel = true;
  int ivar = 0;

  // Preliminary checks
  if (_isKrigingRequested())
  {
    if (_data == nullptr)
    {
      messerr("For this calculation option, you must define some Data");
      return 1;
    }
    if (_data->getLocNumber(ELoc::Z) != 1)
    {
      messerr("The Input dbin must contain ONE variable (Z locator)");
      return 1;
    }
  }
  if (_isSimulationRequested())
  {
    if (nbsimu < 0)
    {
      messerr("For this option, you must define a positive number of simulations");
      return 1;
    }
  }

  if (_isSimulationRequested())
    law_set_random_seed(seed);

  if (_isKrigingRequested())
    _precisionsKrig->makeReady();
  if (_isSimulationRequested())
    _precisionsSimu->makeReady();

  // Preliminary tasks
  if (_data != nullptr)
  {
    dataVect = _data->getColumnByLocator(ELoc::Z,ivar,useSel);
    // Suppress any TEST value and center by the drift
    dataVect = VH::suppressTest(dataVect);
    _centerByDrift(dataVect,ivar,useSel);
  }

  // Create the output vectors in the output Db
  int ncols = 0;
  if (_calcul == ESPDECalcMode::KRIGING)
    ncols = 1;
  else if (_calcul == ESPDECalcMode::KRIGVAR)
    ncols = 2;
  else
    ncols = nbsimu;
  if (ncols <= 0)
  {
    messerr("The number of output attributes should be positive");
    return 1;
  }
  int iptr = dbout->addColumnsByConstant(ncols);

  // Dispatch

  VectorDouble result(dbout->getSampleNumber(true));

  if (_calcul == ESPDECalcMode::KRIGING)
  {
    result.fill(0.);
    _workingData = _workingDataInit;
    _computeKriging();
    for (int icov = 0, ncov = (int) _meshingKrig.size(); icov < ncov; icov++)
      _projecLocal(dbout, _meshingKrig[icov], _workingKrig[icov], result);
    _addDrift(dbout, result);
    dbout->setColumnByUID(result, iptr, useSel);
    namconv.setNamesAndLocators(_data, VectorString(), ELoc::Z, 1, dbout, iptr,
                                "estim", ncols);
  }

  if (_calcul == ESPDECalcMode::KRIGVAR)
  {
    // Estimation by Kriging
    result.fill(0.);
    _workingData = _workingDataInit;
    _computeKriging();
    for (int icov = 0, ncov = (int) _meshingKrig.size(); icov < ncov; icov++)
      _projecLocal(dbout, _meshingKrig[icov], _workingKrig[icov], result);
    _addDrift(dbout, result);
    dbout->setColumnByUID(result, iptr, useSel);
    namconv.setNamesAndLocators(_data, VectorString(), ELoc::Z, 1, dbout, iptr,
                                "estim", 1);


    // Standard Deviation using Monte-Carlo simulations
    VectorDouble temp_mean(dbout->getSampleNumber(true), 0.);
    VectorDouble temp_mean2(dbout->getSampleNumber(true), 0.);

    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      result.fill(0.);
      _computeSimuCond();
      for (int icov = 0, ncov = (int) _meshingSimu.size(); icov < ncov; icov++)
      {
        _projecLocal(dbout, _meshingSimu[icov], _workingSimu[icov], result);
        _projecLocal(dbout, _meshingKrig[icov], _workingKrig[icov], result);
      }
      _addNuggetOnResult(result);
      _addDrift(dbout, result);

      VH::addInPlace(temp_mean, result);
      VH::addSquareInPlace(temp_mean2, result);
    }
    VH::mean1AndMean2ToStdev(temp_mean, temp_mean2, result, nbsimu);
    dbout->setColumnByUID(result, iptr + 1, useSel);
    namconv.setNamesAndLocators(_data, VectorString(), ELoc::Z, 1, dbout, iptr+1,
                                "stdev", 1);
  }

  if (_calcul == ESPDECalcMode::SIMUNONCOND)
  {
    for(int isimu = 0; isimu < nbsimu; isimu++)
    {
      result.fill(0.);
      _computeSimuNonCond();
      for(int icov = 0, ncov = (int) _meshingSimu.size() ; icov < ncov; icov++)
        _projecLocal(dbout, _meshingSimu[icov], _workingSimu[icov], result);
      _addNuggetOnResult(result);
      _addDrift(dbout, result);
      dbout->setColumnByUID(result, iptr + isimu, useSel);
    }
    namconv.setNamesAndLocators(dbout, iptr, "", ncols);
  }

  if (_calcul == ESPDECalcMode::SIMUCOND)
  {
    for(int isimu = 0; isimu < nbsimu; isimu++)
    {
      result.fill(0.);
      _computeSimuCond();
      for(int icov = 0, ncov = (int) _meshingSimu.size(); icov < ncov; icov++)
      {
        _projecLocal(dbout, _meshingSimu[icov], _workingSimu[icov], result);
        _projecLocal(dbout, _meshingKrig[icov], _workingKrig[icov], result);
       }
      _addNuggetOnResult(result);
      _addDrift(dbout, result);
      dbout->setColumnByUID(result, iptr + isimu, useSel);
    }
    namconv.setNamesAndLocators(_data, VectorString(), ELoc::Z, 1, dbout, iptr,
                                "", ncols);
  }
  return iptr;
}

void SPDE::_projecLocal(Db* dbout,
                        const AMesh* meshing,
                        VectorDouble& working,
                        VectorDouble& result)
{
  VectorDouble temp_out(dbout->getSampleNumber(true));
  ProjMatrix proj(dbout,meshing);
  proj.mesh2point(working,temp_out);
  VH::addInPlace(result,temp_out);
}

void SPDE::_addNuggetOnResult(VectorDouble &result) const
{
  if (_nugget <= 0) return;
  for (int iech = 0, nech = (int) result.size(); iech < nech; iech++)
    result[iech] += law_gaussian(0., sqrt(_nugget));
}

bool SPDE::_isSimulationRequested() const
{
  return _calcul == ESPDECalcMode::KRIGVAR
      || _calcul == ESPDECalcMode::SIMUCOND
      || _calcul == ESPDECalcMode::SIMUNONCOND;
}

bool SPDE::_isKrigingRequested() const
{
  return _calcul == ESPDECalcMode::SIMUCOND
      || _calcul == ESPDECalcMode::KRIGING
      || _calcul == ESPDECalcMode::KRIGVAR;
}

double SPDE::computeLogDet(int nbsimu,int seed) const
{
  if (_precisionsKrig == nullptr)
  {
    messerr("The member '_precisionsKrig' must have been calculated beforehand");
    return TEST;
  }

  return _precisionsKrig->computeTotalLogDet(nbsimu,seed);
}

double SPDE::computeQuad() const
{
  if (_data == nullptr)
  {
    messerr("The 'data' must have been spcified beforehand");
    return TEST;
  }
  if (_precisionsKrig == nullptr)
  {
    messerr("The member '_precisionsKrig' must have been calculated beforehand");
    return TEST;
  }

  int ivar = 0;
  bool useSel = true;
  VectorDouble dataVect = _data->getColumnByLocator(ELoc::Z,ivar,useSel);
  _centerByDrift(dataVect,ivar,useSel);
  return _precisionsKrig->computeQuadratic(_workingData);
}

double SPDE::_computeLogLike(int nbsimu, int seed) const
{
  if (_precisionsKrig == nullptr)
  {
    messerr("The member '_precisionsKrig' must have been calculated beforehand");
    return TEST;
  }

  if (!_isCoeffsComputed)
  {
    _computeDriftCoeffs();
  }
  return - 0.5 * (computeLogDet(nbsimu,seed) + computeQuad());
}

/**
 * Calculate the Log-Likelihood profiling the Drift parameters
 */
double SPDE::computeLogLike(int nbsimu, int seed) const
{
  VectorDouble dataVect;
  bool useSel = true;
  int ivar = 0;

  // Preliminary checks
  if (_isKrigingRequested())
  {
    if (_data == nullptr)
    {
      messerr("For this calculation option, you must define some Data");
      return 1;
    }
    if (_data->getLocNumber(ELoc::Z) != 1)
    {
      messerr("The Input dbin must contain ONE variable (Z locator)");
      return 1;
    }
  }

  if (_isKrigingRequested())
    _precisionsKrig->makeReady();

  // Preliminary tasks
  if (_data != nullptr)
  {
    dataVect = _data->getColumnByLocator(ELoc::Z,ivar,useSel);
    // Suppress any TEST value and center by the drift
    dataVect = VH::suppressTest(dataVect);
    _centerByDrift(dataVect,ivar,useSel);
  }

  // Dispatch

  _workingData = _workingDataInit;
  _computeKriging();

  // we assume that covariance parameters have changed when using this function:
  // so driftCoeffs have to be recomputed
  _isCoeffsComputed = false;

  return _computeLogLike(nbsimu,seed);
}

void SPDE::_computeDriftCoeffs() const
{
  if (!_isCoeffsComputed)
   {
    if (_requireCoeffs)
    {
      _precisionsKrig->makeReady();
      _driftCoeffs = _precisionsKrig->computeCoeffs(_data->getColumnByLocator(ELoc::Z,0,true),
                                                    _driftTab);
    }
    _isCoeffsComputed = true;
  }
}

void SPDE::setDriftCoeffs(const VectorDouble& coeffs)
{
  _driftCoeffs  = coeffs;
  _isCoeffsComputed = true;
}

VectorDouble SPDE::getCoeffs()
{
  _computeDriftCoeffs();
  return _driftCoeffs;
}

/**
 * Perform the estimation by KRIGING under the SPDE framework
 *
 * @param dbin Input Db (must contain the variable to be estimated)
 * @param dbout Output Db where the estimation must be performed
 * @param model Model definition
 * @param flag_est True for the estimation
 * @param flag_std True for the standard deviation of estimation error
 * @param mesh Mesh description (optional)
 * @param useCholesky Define the choice regarding Cholesky
 * @param params Set of parameters
 * @param nbMC Number of Monte-Carlo simulations used for variance calculation
 * @param seed Seed used for the Random Number generator
 * @param verbose Verbose flag
 * @param showStats Show statistics for Linear Operations
 * @param namconv Naming convention
 * @return Error return code
 *
 * @remarks You can provide an already existing mesh. Otherwise an optimal mesh will be created
 * @remarks internally: one per structure constituting the Model for Kriging.
 * @remarks Each mesh is created using the Turbo Meshing facility... based on an internal grid.
 * @remarks This internal grid is rotated according to the rotation of the structure. Its mesh size
 * @remarks is derived from the range (per direction) by dividing it by the refinement factor.
 *
 * @remarks Note that switching 'flag_std' to ON implies that 'flag_est' is ON.
 */
int krigingSPDE(Db *dbin,
                Db *dbout,
                Model *model,
                bool flag_est,
                bool flag_std,
                const AMesh *mesh,
                int useCholesky,
                const SPDEParam& params,
                int nbMC,
                int seed,
                bool verbose,
                bool showStats,
                const NamingConvention &namconv)
{
  const ESPDECalcMode mode = (flag_std) ?
      ESPDECalcMode::KRIGVAR : ESPDECalcMode::KRIGING;
  SPDE spde(model, dbout, dbin, mode, mesh, useCholesky, params, verbose,
            showStats);
  return spde.compute(dbout, nbMC, seed, namconv);
}

/**
 * Perform simulations under the SPDE framework
 *
 * @param dbin Input Db. If defined, the simulations are conditional; non conditional otherwise
 * @param dbout Output Db where the simulations must be performed
 * @param model Model definition
 * @param nbsimu Number of simulations
 * @param mesh Mesh description (optional)
 * @param useCholesky Define the choice regarding Cholesky
 * @param params Set of parametes
 * @param seed Seed used for the Random Number generator
 * @param verbose Verbose flag
 * @param showStats Show statistics for Linear Operations
 * @param namconv Naming convention
 * @return Error return code
 *
 * @remarks You can provide an already existing mesh. Otherwise an optimal mesh will be created
 * @remarks internally: one per structure constituting the Model for Kriging and one for Simulating
 * @remarks Each mesh is created using the Turbo Meshing facility... based on an internal grid.
 * @remarks This internal grid is rotated according to the rotation of the structure. Its mesh size
 * @remarks is derived from the range (per direction) by dividing it by the refinement factor.
 */
int simulateSPDE(Db *dbin,
                 Db *dbout,
                 Model *model,
                 int nbsimu,
                 const AMesh *mesh,
                 int useCholesky,
                 const SPDEParam& params,
                 int seed,
                 bool verbose,
                 bool showStats,
                 const NamingConvention &namconv)
{
  const ESPDECalcMode mode = (dbin == nullptr) ?
      ESPDECalcMode::SIMUNONCOND : ESPDECalcMode::SIMUCOND;
  SPDE spde(model, dbout, dbin, mode, mesh, useCholesky, params, verbose,
            showStats);
  return spde.compute(dbout, nbsimu, seed, namconv);
}

double logLikelihoodSPDE(Db *dbin,
                         Db *dbout,
                         Model *model,
                         const AMesh *mesh,
                         int useCholesky,
                         int nbsimu,
                         int seed,
                         const SPDEParam& params,
                         bool verbose)
{
  SPDE spde(model, dbout, dbin, ESPDECalcMode::KRIGING, mesh, useCholesky,
            params, verbose, false);
  return spde.computeLogLike(nbsimu, seed);
}

static int _loadPositions(int iech,
                          const VectorVectorInt &index1,
                          const VectorInt &cumul,
                          VectorInt& positions,
                          VectorInt& identity,
                          int *rank_arg)
{
  int nvar = (int) cumul.size();
  int ndef = 0;
  int rank = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    rank = 2 * rank;
    int ipos = VH::whereElement(index1[ivar], iech);
    if (ipos < 0)
      positions[ivar] = -1;
    else
    {
      positions[ivar] = ipos + cumul[ivar];
      identity[ndef] = ivar;
      ndef++;
      rank += 1;
    }
  }
  *rank_arg = rank;
  return ndef;
}

static void _addVerrConstant(MatrixSquareSymmetric& sills, const VectorDouble& verrDef)
{
  int nverr = (int) verrDef.size();
  if (nverr > 0)
  {
    for (int iverr = 0; iverr < nverr; iverr++)
      sills.updValue(iverr, iverr, EOperator::ADD, verrDef[iverr]);
  }
}

static void _checkMinNugget(MatrixSquareSymmetric& sills, const VectorDouble& minNug)
{
  int nvar = (int) minNug.size();

  // Check that the diagonal of the Sill matrix is large enough
  for (int ivar = 0; ivar < nvar; ivar++)
    sills.setValue(ivar,ivar, MAX(sills.getValue(ivar,ivar), minNug[ivar]));
}

static MatrixSquareSymmetric _getSillGlobalMatrix(Model *model, int icovNug)
{
  int nvar = model->getVariableNumber();

  // Elaborate the matrix of sills for the Nugget Effect component
   MatrixSquareSymmetric sills(nvar);
   if (icovNug >= 0) sills = model->getSillValues(icovNug);
   return sills;
}

static MatrixSquareSymmetric _buildSillPartialMatrix(const MatrixSquareSymmetric &sillsRef,
                                                     int nvar,
                                                     int ndef,
                                                     const VectorInt &identity)
{
  MatrixSquareSymmetric sills;
  if (ndef == nvar)
    sills = sillsRef;
  else
  {
    sills = MatrixSquareSymmetric(ndef);
    for (int idef = 0; idef < ndef; idef++)
      for (int jdef = 0; jdef <= idef; jdef++)
        sills.setValue(idef, jdef, sillsRef.getValue(identity[idef], identity[jdef]));
  }
  return sills;
}

/**
 * Build the inverse of the Nugget Effect matrix
 * It is established for:
 * - the number of variables defined in 'dbin' (and in 'Model')
 * - the active samples of 'dbin'
 * - the samples where Z-variable (and possibly V-variable) is defined
 *
 * @param db Input Db structure
 * @param model Input Model
 * @param params SPDEParam structure
 */
MatrixSparse* buildInvNugget(Db *db, Model *model, const SPDEParam& params)
{
  MatrixSparse* mat = nullptr;
  if (db == nullptr) return mat;
  int nech = db->getSampleNumber();
  if (model == nullptr) return mat;
  int nvar = db->getLocNumber(ELoc::Z);
  if (nvar != model->getVariableNumber())
  {
    messerr("'db' and 'model' should have the same number of variables");
    return mat;
  }
  VectorInt ivars = VH::sequence(nvar);

  // Get the minimum value for diagonal terms
  double eps = params.getEpsNugget();
  VectorDouble minNug(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    minNug[ivar] = eps * model->getTotalSill(ivar, ivar);

  // Play the non-stationarity (if needed)
  ANoStat *nostat = model->getNoStatModify();
  bool flag_nostat_sill = (nostat != nullptr && nostat->isDefinedByType(EConsElem::SILL));
  if (flag_nostat_sill)
  {
    if (nostat->manageInfo(1, db, db) != 0) return mat;
  }

  // Create the sets of Vector of valid sample indices per variable (not masked and defined)
  VectorVectorInt index1 = db->getMultipleRanksActive(ivars);
  // 'cumul' counts the number of valid positions for all variables before 'ivar'
  VectorInt cumul(nvar, 0);
  int number = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    cumul[ivar] = number;
    number += (int) index1[ivar].size();
  }

  // Check the various possibilities
  // - flag_verr: True if Variance of Measurement Error variable is defined
  // - flag_isotropic: True in Isotopic case
  // - flag_uniqueVerr: True if the Variance of Measurement Error is constant per variable
  // - flag_nostat: True is some non-stationarity is defined
  int nverr = db->getLocNumber(ELoc::V);
  bool flag_verr = (nverr > 0);
  bool flag_isotopic = true;
  for (int ivar = 1; ivar < nvar && flag_isotopic; ivar++)
    if (! VH::isSame(index1[ivar], index1[0])) flag_isotopic = false;
  bool flag_uniqueVerr = true;
  VectorDouble verrDef(nverr, 0.);
  if (flag_verr)
  {
    for (int iverr = 0; iverr < nverr && flag_uniqueVerr; iverr++)
    {
      VectorDouble verr = db->getColumnByLocator(ELoc::V, iverr);
      if ((int) VH::unique(verr).size() > 1) flag_uniqueVerr = false;
      verrDef[iverr] = verr[0];
    }
  }
  bool flag_constant = (! flag_nostat_sill && (! flag_verr || flag_uniqueVerr));

  // Elaborate the Sill matrix for the Nugget Effect component
  int icovNug = model->getRankNugget();
  MatrixSquareSymmetric sillsRef = _getSillGlobalMatrix(model, icovNug);
  int count = (int) pow(2, nvar);
  std::vector<MatrixSquareSymmetric> sillsInv(count);

  // Pre-calculate the inverse of the sill matrix (if constant)

  if (flag_constant)
  {
    // In case of (Unique) Variance of measurement error, patch sill matrix
    if (flag_verr) _addVerrConstant(sillsRef, verrDef);

    // Check that the diagonal of the Sill matrix is large enough
    _checkMinNugget(sillsRef, minNug);
  }

  // Constitute the triplet
  NF_Triplet NF_T;

  // Loop on the samples
  int rank;
  int ndef = nvar;
  VectorInt position(nvar);
  VectorInt identity(nvar);
  for (int iech = 0; iech < nech; iech++)
  {
    if (! db->isActive(iech)) continue;

    // Count the number of variables for which current sample is valid
    ndef = _loadPositions(iech, index1, cumul, position, identity, &rank);
    if (ndef <= 0) continue;

    // If all samples are defined, in the stationary case, use the inverted sill matrix
    if (flag_constant)
    {
      if (sillsInv[rank].empty())
      {
        sillsInv[rank] = _buildSillPartialMatrix(sillsRef, nvar, ndef, identity);
        if (sillsInv[rank].invert() != 0) return mat;
      }
      for (int idef = 0; idef < ndef; idef++)
        for (int jdef = 0; jdef < ndef; jdef++)
          NF_T.add(position[identity[idef]], position[identity[jdef]],
                   sillsInv[rank].getValue(idef,jdef));
    }
    else
    {
      // Update due to non-stationarity (optional)
      if (flag_nostat_sill)
      {
        model->updateCovByPoints(1, iech, 1, iech);
        sillsRef = _getSillGlobalMatrix(model, icovNug);
      }

      // Establish a local matrix
      MatrixSquareSymmetric local(ndef);
      for (int idef = 0; idef < ndef; idef++)
        for (int jdef = 0; jdef <= idef; jdef++)
        {
          // Load the sill value of the Nugget Effect component
          double value = sillsRef.getValue(identity[idef], identity[jdef]);

          // Patch the diagonal term of the local matrix
          if (idef == jdef)
          {
            // Add the Variance of measurement error (optional)
            if (flag_verr && idef < nverr)
              value += db->getFromLocator(ELoc::V, iech, identity[idef]);

            // Check the minimum values over the diagonal
            value = MAX(value, MAX(local.getValue(idef, idef), minNug[idef]));
          }

          local.setValue(idef, jdef, value);
        }
      if (local.invert() != 0) return mat;

      for (int idef = 0; idef < ndef; idef++)
         for (int jdef = 0; jdef < ndef; jdef++)
           NF_T.add(position[identity[idef]],position[identity[jdef]],
                    local.getValue(idef,jdef));
    }
  }

  // Convert from triplet to sparse matrix
  mat = MatrixSparse::createFromTriplet(NF_T);

  // Free the non-stationary specific allocation
  if (model->isNoStat())
    (void) nostat->manageInfo(-1, db, db);

  return mat;
}
