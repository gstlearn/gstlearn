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
#include "Covariances/CovAniso.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/AException.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/NamingConvention.hpp"
#include "Model/Model.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

#include <iostream>
#include <math.h>

/**
 * The class constructor with the following arguments:
 *
 * @param model  This compulsory argument is a LMC of Matern's (or Markov?) basic structures with possibly a nugget effect
 * @param domain The Db defining the space dimension and spatial limits where the SPDE model is defined
 * @param data   The Db containing the data for conditioning (optional)
 * @param calcul Option from ESPDECalcMode
 * @param meshUser The mesh for the discretization of the domain
 * @param verbose  Verbose flag
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
           bool verbose)
    : _data(nullptr),
      _calcul(),
      _refineK(11),
      _refineS(18),
      _border(8),
      _precisionsKrig(nullptr),
      _precisionsSimu(nullptr),
      _pilePrecisions(),
      _pileProjMatrix(),
      _meshingSimu(),
      _meshingKrig(),
      _driftCoeffs(),
      _model(nullptr),
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
      _nIterMax(1000),
      _eps(EPSILON8),
      _epsNugget(1.e-2)
{
  (void) _init(model, domain, data, calcul, meshUser, verbose);
}

SPDE::~SPDE()
{
  //_purge();
}

void SPDE::_purge()
{
  delete _precisionsKrig;
  delete _precisionsSimu;

  for(auto &e : _pilePrecisions)
  {
    delete e;
  }
  _pilePrecisions.clear();

  for(auto &e : _pileProjMatrix)
  {
    delete e;
  }
  _pileProjMatrix.clear();
  for(auto &e : _projOnDbOut)
  {
    delete e;
  }

  if (_deleteMesh)
  {
    for(auto &e : _meshingSimu)
    {
      delete e;
    }
    for(auto &e : _meshingKrig)
    {
      delete e;
    }
  }
  _meshingSimu.clear();
  _meshingKrig.clear();
}

SPDE* SPDE::create(Model *model,
                   const Db *domain,
                   const Db *data,
                   const ESPDECalcMode &calcul,
                   const AMesh* meshUser,
                   bool verbose)
{
  return new SPDE(model, domain, data, calcul, meshUser, verbose);
}

bool SPDE::_useCholesky() const
{
  if (_model->getDimensionNumber() == 2)
    return true;
  else
    return false;
}

int SPDE::_init(Model *model,
                const Db *domain,
                const Db *data,
                const ESPDECalcMode &calcul,
                const AMesh *meshUser,
                bool verbose)
{
  //_purge();
  _model  = model;
  _calcul = calcul;
  _data   =  data;
  const ANoStat* nostat = nullptr;
  if (model->isNoStat()) nostat = model->getNoStat();

  // Preliminary check
  if (_performKriging() && _data == nullptr)
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
  if (_performSimulation())
    _precisionsSimu = new PrecisionOpMultiConditional();
  if (_performKriging() || _requireCoeffs)
    _precisionsKrig = new PrecisionOpMultiConditional();

  // Loop on the basic structures
  for(int icov = 0 ; icov < model->getCovaNumber(); icov++)
  {
    const CovAniso* cova = model->getCova(icov);
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

      if (_performSimulation())
      {
        if (meshUser == nullptr)
        {
          mesh = MeshETurbo::createFromCova(*cova, domain, _refineS, _border,
                                            useSel, flagNoStatRot, verbose);
          _deleteMesh = true;
        }

        if (_useCholesky())
          precision = new PrecisionOpCs(mesh, model, icov, verbose);
        else
          precision = new PrecisionOp(mesh, model, icov, verbose);
        _pilePrecisions.push_back(precision);

        proj = new ProjMatrix(_data, mesh, verbose);
        _pileProjMatrix.push_back(proj);

        _meshingSimu.push_back(mesh);

        _precisionsSimu->push_back(precision,proj);
        _precisionsSimu->setVarianceDataVector(varianceData);
        _workingSimu.push_back(VectorDouble(precision->getSize()));
      }

      if (_performKriging() || _requireCoeffs)
      {
        if (meshUser == nullptr)
        {
          mesh = MeshETurbo::createFromCova(*cova, domain, _refineK, _border,
                                            useSel, flagNoStatRot, verbose);
          _deleteMesh = true;
        }
        _meshingKrig.push_back(mesh);

        if (_useCholesky())
          precision = new PrecisionOpCs(mesh, model, icov);
        else
          precision = new PrecisionOp(mesh, model, icov);
        _pilePrecisions.push_back(precision);

        proj = new ProjMatrix(_data,mesh);
        _pileProjMatrix.push_back(proj);

        _precisionsKrig->setNIterMax(_nIterMax);
        _precisionsKrig->setEps(_eps);
        _precisionsKrig->push_back(precision,proj);
        _workingKrig.push_back(VectorDouble(precision->getSize()));
      }
    }
    else
    {
      messerr("SPDE is only implemented for Matérn covariances (BESSEL_K) and Markov (MARKOV)");
      return 1;
    }
  }

  // Evaluation of the variance at data point
  // (nugget + measurement error or minimum proportion of total sill)
  if (_performKriging())
  {
    if (_data->getLocNumber(ELoc::V) > 0)
    {
      varianceData = _data->getColumnByLocator(ELoc::V,0,useSel);
      for (int iech = 0; iech < _data->getSampleNumber(true); iech++)
      {
        double *temp_dat = &varianceData[iech];
        *temp_dat = MAX(*temp_dat + _nugget, _epsNugget * totalSill);
      }
    }
    else
    {
      VH::fill(varianceData, MAX(_nugget, _epsNugget * totalSill),
               _data->getSampleNumber(true));
    }
    _precisionsKrig->setVarianceDataVector(varianceData);

    if (_performSimulation())
    {
      _precisionsSimu->setVarianceDataVector(varianceData);
    }
  }
  return 0;
}

void SPDE::_computeLk() const
{
  VectorVectorDouble rhs = _precisionsKrig->computeRhs(_workingData);
  _precisionsKrig->initLk(rhs,_workingKrig); // Same as evalInverse but with just one iteration
}

void SPDE::_computeKriging() const
{
  VectorVectorDouble rhs = _precisionsKrig->computeRhs(_workingData);
  _precisionsKrig->evalInverse(rhs,_workingKrig);
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
    _workingDataInit = _model->evalDrifts(_data,_driftCoeffs,ivar,useSel);

    for(int iech = 0, nech = (int) _workingDataInit.size(); iech<nech; iech++)
    {
      _workingDataInit[iech] = dataVect[iech] - _workingDataInit[iech];
    }
  }
}

void SPDE::_addDrift(Db* db, VectorDouble &result, int ivar, bool useSel)
{
  if (! _requireCoeffs) return;
  VectorDouble temp_out = _model->evalDrifts(db, _driftCoeffs, ivar, useSel);
  VH::addInPlace(result, temp_out);
}

int SPDE::compute(Db *dbout, int nbsimu, int seed, const NamingConvention &namconv)
{
  VectorDouble dataVect;
  bool useSel = true;
  int ivar = 0;

  // Preliminary checks
  if (_performKriging())
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
  if (_performSimulation() && nbsimu < 0)
  {
    messerr("For this option, you must define a positive number of simulations");
    return 1;
  }

  // Preliminary tasks
  if (_data != nullptr)
  {
    dataVect = _data->getColumnByLocator(ELoc::Z,ivar,useSel);
    _centerByDrift(dataVect,ivar,useSel);
  }

  if (_performSimulation())
    law_set_random_seed(seed);

  // Add the vectors in the output Db
  int ncols = 1;
  if (_performSimulation()) ncols = nbsimu;
  int iptr = dbout->addColumnsByConstant(ncols);

  // Dispatch

  VectorDouble temp_out(dbout->getSampleNumber(true));
  String suffix = "";
  bool useVarName = false;

  if (_calcul == ESPDECalcMode::KRIGING)
  {
    VectorDouble result(dbout->getSampleNumber(true),0.);
    _workingData = _workingDataInit;
    _computeKriging();
    for(int icov = 0, ncov = (int) _meshingKrig.size() ; icov < ncov; icov++)
    {
      ProjMatrix projKriging(dbout,_meshingKrig[icov]);
      projKriging.mesh2point(_workingKrig[icov],temp_out);
      VH::addInPlace(result,temp_out);
    }
    _addDrift(dbout, result);
    dbout->setColumnByUID(result, iptr, useSel);
    useVarName = true;
    suffix = "estim";
  }

  if (_calcul == ESPDECalcMode::SIMUNONCOND)
  {
    for(int isimu = 0; isimu < nbsimu; isimu++)
    {
      VectorDouble result(dbout->getSampleNumber(true),0.);
      _computeSimuNonCond();
      for(int icov = 0, ncov = (int) _meshingSimu.size() ; icov < ncov; icov++)
      {
        ProjMatrix projSimu(dbout,_meshingSimu[icov]);
        projSimu.mesh2point(_workingSimu[icov],temp_out);
        VH::addInPlace(result,temp_out);
      }
      _addNuggetOnResult(result);
      _addDrift(dbout, result);
      dbout->setColumnByUID(result, iptr + isimu, useSel);
    }
  }

  if (_calcul == ESPDECalcMode::SIMUCOND)
  {
    for(int isimu = 0; isimu < nbsimu; isimu++)
    {
      VectorDouble result(dbout->getSampleNumber(true),0.);
      _computeSimuCond();
      for(int icov = 0, ncov = (int) _meshingSimu.size(); icov < ncov; icov++)
      {
        ProjMatrix projSimu(dbout,_meshingSimu[icov]);
        projSimu.mesh2point(_workingSimu[icov],temp_out);
        VH::addInPlace(result,temp_out);
        ProjMatrix projKriging(dbout,_meshingKrig[icov]);
        projKriging.mesh2point(_workingKrig[icov],temp_out);
        VH::addInPlace(result,temp_out);
       }
      _addNuggetOnResult(result);
      _addDrift(dbout, result);
      dbout->setColumnByUID(result, iptr + isimu, useSel);
    }
    useVarName = true;
  }

  if (_calcul == ESPDECalcMode::LIKELIHOOD)
  {
    VectorDouble result(dbout->getSampleNumber(true),0.);
    _computeLk();
    for(int icov = 0, ncov = (int) _meshingKrig.size() ; icov < ncov; icov++)
    {
      ProjMatrix proj(dbout,_meshingKrig[icov]);
      proj.mesh2point(_workingKrig[icov],temp_out);
      VH::addInPlace(result,temp_out);
    }
    _addDrift(dbout, result);
    dbout->setColumnByUID(result, iptr, useSel);
    suffix = "likelihood";
  }

  if (useVarName)
    namconv.setNamesAndLocators(_data, ELoc::Z, 1, dbout, iptr, suffix, ncols);
  else
    namconv.setNamesAndLocators(dbout, iptr, suffix, ncols);
  return iptr;
}

void SPDE::_addNuggetOnResult(VectorDouble &result)
{
  if (_nugget <= 0) return;
  for (int iech = 0, nech = (int) result.size(); iech < nech; iech++)
    result[iech] += law_gaussian(0., sqrt(_nugget));
}

bool SPDE::_performSimulation() const
{
  return _calcul == ESPDECalcMode::SIMUCOND
      || _calcul == ESPDECalcMode::SIMUNONCOND;
}

bool SPDE::_performKriging() const
{
  return _calcul == ESPDECalcMode::SIMUCOND
      || _calcul == ESPDECalcMode::KRIGING
      || _calcul == ESPDECalcMode::LIKELIHOOD;
}

double SPDE::computeLogDet(int nbsimu,int seed) const
{
  double val;
  val = _precisionsKrig->computeTotalLogDet(nbsimu,seed);
  return val;
}

double SPDE::computeQuad() const
{
  if (_data == nullptr) return TEST;
  int ivar = 0;
  bool useSel = true;
  VectorDouble dataVect = _data->getColumnByLocator(ELoc::Z,ivar,useSel);
  _centerByDrift(dataVect,ivar,useSel);
  return _precisionsKrig->computeQuadratic(_workingData);
}

double SPDE::computeLogLike(int nbsimu, int seed) const
{
  if (!_isCoeffsComputed)
  {
    _computeDriftCoeffs();
  }
  return - 0.5 * (computeLogDet(nbsimu,seed) + computeQuad()) ;
}

double SPDE::computeProfiledLogLike(int nbsimu, int seed) const
{
  _isCoeffsComputed = false; // we assume that covariance parameters have changed when using this function
  //  so driftCoeffs have to be recomputed

  return computeLogLike(nbsimu,seed);
}

void SPDE::_computeDriftCoeffs() const
{
  if (!_isCoeffsComputed)
   {
    _isCoeffsComputed = true;
    if (_requireCoeffs)
    {
      _driftCoeffs = _precisionsKrig->computeCoeffs(_data->getColumnByLocator(ELoc::Z,0,true),
                                                    _driftTab);
    }
  }
}

void SPDE::setDriftCoeffs(VectorDouble coeffs)
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
 * @param flag_varz True for the variance of the estimator
 * @param mesh Mesh description (optional)
 * @param refineK Refinement factor for building internal meshing for Kriging
 * @param border Number of nodes used for extending the internal grid
 * @param epsNugget Value for the relaxing nugget effect
 * @param verbose Verbose flag
 * @param namconv Naming convention
 * @return Error return code
 *
 * @remarks You can provide an already existing mesh. Otherwise an optimal mesh will be created
 * @remarks internally: one per structure constituting the Model for Kriging.
 * @remarks Each mesh is created using the Turbo Meshing facility... based on an internal grid.
 * @remarks This internal grid is rotated according to the rotation of the structure. Its mesh size
 * @remarks is derived from the range (per direction) by dividing it by the refinement factor.
 */
int krigingSPDE(Db *dbin,
                Db *dbout,
                Model *model,
                bool flag_est,
                bool flag_std,
                bool flag_varz,
                const AMesh *mesh,
                int refineK,
                int border,
                double epsNugget,
                bool verbose,
                const NamingConvention &namconv)
{
  // Preliminary checks
  if (flag_std || flag_varz)
  {
    messerr("These options have not been implemented yet. Not taken into account");
  }
  SPDE spde(model, dbout, dbin, ESPDECalcMode::KRIGING, mesh, verbose);
  spde.setRefineK(refineK);
  spde.setBorder(border);
  spde.setEpsNugget(epsNugget);

  return spde.compute(dbout, 0, 0, namconv);
}

/**
 * Perform simulations under the SPDE framework
 *
 * @param dbin Input Db. If defined, the simulations are conditional; non conditional otherwise
 * @param dbout Output Db where the simulations must be performed
 * @param model Model definition
 * @param nbsimu Number of simulations
 * @param mesh Mesh description (optional)
 * @param refineK Refinement factor for building internal meshing for Kriging
 * @param refineS Refinement factor for building internal meshing for Simulations
 * @param border Number of nodes used for extending the internal grid
 * @param seed Seed used for the Random Number generator
 * @param epsNugget Value for the relaxing nugget effect
 * @param verbose Verbose flag
 * @param namconv Naming convention
 * @return Error return code
 *
 * @remarks You can provide an already existing mesh. Otherwise an optimal mesh will be created
 * @remarks internally: one per structure constituting the Model for Kriging and one for Simulating
 * @remarks Each mesh is created using the Turbo Meshing facility... based on an internal grid.
 * @remarks This internal grid is rotated according to the rotation of the structure. Its mesh size
 * @remarks is derived from the range (per direction) by dividing it by the refinement factor.
 */
GSTLEARN_EXPORT int simulateSPDE(Db *dbin,
                                 Db *dbout,
                                 Model *model,
                                 int nbsimu,
                                 const AMesh *mesh,
                                 int refineK,
                                 int refineS,
                                 int border,
                                 int seed,
                                 double epsNugget,
                                 bool verbose,
                                 const NamingConvention &namconv)
{
  const ESPDECalcMode mode = (dbin == nullptr) ? ESPDECalcMode::SIMUNONCOND : ESPDECalcMode::SIMUCOND;
  SPDE spde(model, dbout, dbin, mode, mesh, verbose);
  spde.setRefineK(refineK);
  spde.setRefineS(refineS);
  spde.setBorder(border);
  spde.setEpsNugget(epsNugget);

  return spde.compute(dbout, nbsimu, seed, namconv);
}

