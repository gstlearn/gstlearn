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

SPDE::SPDE()
    : _data(nullptr),
      _calcul(),
      _precisionsKriging(nullptr),
      _precisionsSimu(nullptr),
      _pilePrecisions(),
      _pileProjMatrix(),
      _simuMeshing(),
      _krigingMeshing(),
      _driftCoeffs(),
      _model(nullptr),
      _workKriging(),
      _workingSimu(),
      _workingData(),
      _workingData2(),
      _projOnDbOut(),
      _adressesICov(),
      _nugget(0.),
      _driftTab(),
      _requireCoeffs(false),
      _isCoeffsComputed(false),
      _deleteMesh(false),
      _nIterMax(1000),
      _eps(EPSILON8)
{
}

SPDE::SPDE(Model* model,
           const DbGrid* field,
           const Db* dat,
           const ESPDECalcMode& calc)
    : _data(nullptr),
      _calcul(),
      _precisionsKriging(nullptr),
      _precisionsSimu(nullptr),
      _pilePrecisions(),
      _pileProjMatrix(),
      _simuMeshing(),
      _krigingMeshing(),
      _driftCoeffs(),
      _model(nullptr),
      _workKriging(),
      _workingSimu(),
      _workingData(),
      _workingData2(),
      _projOnDbOut(),
      _adressesICov(),
      _nugget(0.),
      _driftTab(),
      _requireCoeffs(false),
      _isCoeffsComputed(false),
      _deleteMesh(false),
      _nIterMax(1000),
      _eps(EPSILON8)
{
  init(model, field, dat, calc);
}

SPDE::~SPDE()
{
  //_purge();
}

void SPDE::_purge()
{
  delete _precisionsKriging;
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
    for(auto &e : _simuMeshing)
    {
      delete e;
    }
    for(auto &e : _krigingMeshing)
    {
      delete e;
    }
  }
  _simuMeshing.clear();
  _krigingMeshing.clear();
}

SPDE* SPDE::create(Model *model,
                   const DbGrid *field,
                   const Db *data,
                   const ESPDECalcMode &calc)
{
  return new SPDE(model, field, data, calc);
}

void SPDE::init(Model* model,
                const DbGrid* field,
                const Db* data,
                const ESPDECalcMode& calc,
                const AMesh* meshUser,
                bool verbose)
{
  //_purge();
  _model  = model;
  _calcul = calc;
  _data   =  data;
  bool useSel = true;
  VectorDouble varianceData;
  double totalSill = 0.;
  PrecisionOp* precision;
  const AMesh* mesh = meshUser;
  ProjMatrix* proj;
  if (_data != nullptr) _driftTab = _model->getDrifts(_data, useSel);
  _requireCoeffs = _driftTab.size()>0 && _data != nullptr;

  // Loop on the basic structures
  for(int icov = 0 ; icov < model->getCovaNumber(); icov++)
  {
    const CovAniso* cova = model->getCova(icov);

    if (cova->getType() == ECov::NUGGET)
    {
      _nugget = cova->getSill(0,0);
    }
    else if (cova->getType() == ECov::BESSEL_K || cova->getType() == ECov::MARKOV)
    {
      totalSill += cova->getSill(0, 0);
      if (_calculSimu())
      {
        if (meshUser == nullptr)
        {
          mesh = new MeshETurbo();
          ((MeshETurbo*)mesh)->initFromCova(*cova,field,18,5,useSel,verbose);
          _deleteMesh = true;
        }

        if (model->getDimensionNumber() == 2)
        {
          precision = new PrecisionOpCs(mesh, model, icov, verbose);
        }
        else
        {
          precision = new PrecisionOp(mesh, model, icov, verbose);
        }

        _pilePrecisions.push_back(precision);

        proj = new ProjMatrix(_data, mesh, verbose);

        _pileProjMatrix.push_back(proj);

        _simuMeshing.push_back(mesh);

        _precisionsSimu = new PrecisionOpMultiConditional();
        _precisionsSimu->push_back(precision,proj);
        _precisionsSimu->setVarianceDataVector(varianceData);

        _workingSimu.push_back(VectorDouble(precision->getSize()));
      }
      if (_calculKriging() || _requireCoeffs)
      {
        _workingData2.resize(_data->getSampleNumber(useSel));
        if (meshUser == nullptr)
        {
          mesh = new MeshETurbo();
          ((MeshETurbo*)mesh)->initFromCova(*cova,field,11,5,useSel,verbose);
          _deleteMesh = true;
        }
        _krigingMeshing.push_back(mesh);

        if (model->getDimensionNumber() == 2)
        {
          precision = new PrecisionOpCs(mesh, model, icov);
        }
        else
        {
          precision = new PrecisionOp(mesh, model, icov);
        }

        _pilePrecisions.push_back(precision);

        proj = new ProjMatrix(_data,mesh);
        _pileProjMatrix.push_back(proj);

        _precisionsKriging = new PrecisionOpMultiConditional();
        _precisionsKriging->setNIterMax(_nIterMax);
        _precisionsKriging->setEps(_eps);
        _precisionsKriging->push_back(precision,proj);
        _workKriging.push_back(VectorDouble(precision->getSize()));
      }
    }
    else
    {
      my_throw("SPDE is only implemented for Matérn covariances (BESSEL_K) and Markov (MARKOV)");
    }
  }

  // Evaluation of the variance at data point
  // (nugget + measurement error or minimum proportion of total sill)
  if (_calculKriging())
  {
    if (_data->getLocNumber(ELoc::V) > 0)
    {
      varianceData = _data->getColumnByLocator(ELoc::V,0,useSel);
      for (int iech = 0; iech < _data->getSampleNumber(true); iech++)
      {
        double *temp = &varianceData[iech];
        *temp = MAX(*temp+_nugget,0.01 * totalSill);
      }
    }
    else
    {
      VH::fill(varianceData, MAX(_nugget, 0.01 * totalSill),
               _data->getSampleNumber(true));
    }
    _precisionsKriging->setVarianceDataVector(varianceData);

    if (_calculSimu())
    {
      _precisionsSimu->setVarianceDataVector(varianceData);
    }
  }
}

void SPDE::computeLk() const
{
  VectorVectorDouble rhs = _precisionsKriging->computeRhs(_workingData);
  _precisionsKriging->initLk(rhs,_workKriging); // Same as evalInverse but with just one iteration
}

void SPDE::computeKriging() const
{
  VectorVectorDouble rhs = _precisionsKriging->computeRhs(_workingData);
  _precisionsKriging->evalInverse(rhs,_workKriging);
}

void SPDE::computeSimuNonCond(int nbsimus, int seed) const
{
  law_set_random_seed(seed);
  VectorDouble gauss;
  VectorDouble resultSimu;
  for(int isim = 0; isim < nbsimus; isim++)
  {
    for(int icov = 0; icov < (int)_simuMeshing.size();icov++)
    {
      gauss = VH::simulateGaussian(_simuMeshing[icov]->getNApices());
      _precisionsSimu->simulateOnMeshing(gauss,_workingSimu);
    }
  }
}

void SPDE::computeSimuCond(int nbsimus, int seed) const
{
  computeSimuNonCond(nbsimus,seed);
  VectorDouble temp(_data->getSampleNumber(true));
  _precisionsSimu->simulateOnDataPointFromMeshings(_workingSimu,temp);
  VH::multiplyConstant(temp,-1.);
  VH::addInPlace(_workingData,temp);
  computeKriging();
}

void SPDE::centerByDrift(const VectorDouble& dataVect,int ivar,bool useSel) const
{
  _computeDriftCoeffs();

  if (_driftCoeffs.empty())
  {
    if (_workingData.empty())
    {
      _workingData.resize(dataVect.size());
    }

    for(int iech = 0; iech<(int)_workingData.size();iech++)
    {
      _workingData[iech] = dataVect[iech];
    }
  }
  else
  {
    _workingData = _model->evalDrifts(_data,_driftCoeffs,ivar,useSel);

    for(int iech = 0; iech<(int)_workingData.size();iech++)
    {
      _workingData[iech] = dataVect[iech] - _workingData[iech];
    }
  }
}

void SPDE::compute(int nbsimus, int seed)
{
  VectorDouble dataVect;
  bool useSel = true;
  int ivar = 0;

  if (_data != nullptr)
  {
    dataVect = _data->getColumnByLocator(ELoc::Z,ivar,useSel);
    centerByDrift(dataVect,ivar,useSel);
  }

  if (_calcul == ESPDECalcMode::KRIGING)
  {
    computeKriging();
  }

  if (_calcul == ESPDECalcMode::SIMUNONCOND)
  {
    computeSimuNonCond(nbsimus,seed);
  }

  if (_calcul == ESPDECalcMode::SIMUCOND)
  {
    computeSimuCond(nbsimus,seed);
  }

  if (_calcul == ESPDECalcMode::LIKELIHOOD)
  {
    computeLk();
  }
}

MeshETurbo* SPDE::_createMeshing(const CovAniso & cova,
                                const DbGrid& field,
                                double discr,
                                double ext)
{
  VectorDouble extendMin,extendMax;
  bool useSel = true;
  int dim = cova.getNDim();
  for(int idim = 0;idim<(int)dim;idim++)
  {
      auto limits = field.getExtrema(idim,useSel);
      extendMin.push_back(limits[0]);
      extendMax.push_back(limits[1]);
  }
  VectorDouble cellSize(dim,VH::minimum(cova.getRanges())/discr);

  VectorInt nx;
  VectorDouble dx;
  VectorDouble x0;
  double delta;
  for(int idim=0;idim<dim;idim++)
  {
    delta = extendMax[idim]-extendMin[idim];
    nx.push_back((int)(( delta * ( 1 + 2 * ext ) ) / cellSize[idim] ));
    x0.push_back(field.getX0(idim)- delta * ext);
  }
  MeshETurbo* mesh = new MeshETurbo(nx,cellSize,x0,field.getRotMat());
  return mesh;
}

bool SPDE::_calculSimu() const
{
  return _calcul == ESPDECalcMode::SIMUCOND
      || _calcul == ESPDECalcMode::SIMUNONCOND;
}

bool SPDE::_calculKriging() const
{
  return ((_calcul == ESPDECalcMode::SIMUCOND
        || _calcul == ESPDECalcMode::KRIGING
        || _calcul == ESPDECalcMode::LIKELIHOOD)
        && _data != nullptr);
}

int SPDE::query(Db* db, const NamingConvention& namconv) const
{
  int ivar = 0;
  bool useSel = true;
  VectorDouble temp(db->getSampleNumber(true));
  VectorDouble result(db->getSampleNumber(true),0.);
  String suffix;
  if(_calcul == ESPDECalcMode::KRIGING || _calcul == ESPDECalcMode::LIKELIHOOD)
  {
    for(int i = 0 ; i< (int)_krigingMeshing.size(); i++)
    {
      ProjMatrix proj(db,_krigingMeshing[i]);
      proj.mesh2point(_workKriging[i],temp);
      VH::addInPlace(result,temp);
    }
    suffix = "kriging";
  }
  else if(_calcul == ESPDECalcMode::SIMUNONCOND)
  {
    for(int i = 0 ; i< (int)_simuMeshing.size(); i++)
    {
      ProjMatrix proj(db,_simuMeshing[i]);
      proj.mesh2point(_workingSimu[i],temp);
      VH::addInPlace(result,temp);
    }
    //TODO check variance
    for(int iech = 0 ; iech< (int)result.size(); iech++)
    {
      result[iech] += law_gaussian(0.,sqrt(_nugget));
    }
    suffix = "simu";
  }
  else if(_calcul == ESPDECalcMode::SIMUCOND)
  {
    for(int i = 0 ; i< (int)_simuMeshing.size(); i++)
    {
      ProjMatrix projSimu(db,_simuMeshing[i]);
      projSimu.mesh2point(_workingSimu[i],temp);
      VH::addInPlace(result,temp);
      ProjMatrix projKriging(db,_krigingMeshing[i]);
      projKriging.mesh2point(_workKriging[i],temp);
      VH::addInPlace(result,temp);
      // TODO add nugget
     }
    suffix = "condSimu";
  }

  if (_requireCoeffs)
  {
    temp = _model->evalDrifts(db,_driftCoeffs,ivar,useSel);
    VH::addInPlace(result,temp);
  }
  int iptr = db->addColumns(result,"SPDE",ELoc::Z,0,useSel);
  namconv.setNamesAndLocators(_data,ELoc::Z,1,db,iptr,suffix);
  return iptr;
}

double SPDE::computeLogDet(int nbsimus,int seed) const
{
  double val;
  val = _precisionsKriging->computeTotalLogDet(nbsimus,seed);
  return val;
}

double SPDE::computeQuad() const
{
  if (_data == nullptr) return TEST;
  int ivar = 0;
  bool useSel = true;
  VectorDouble dataVect = _data->getColumnByLocator(ELoc::Z,ivar,useSel);
  centerByDrift(dataVect,ivar,useSel);
  return _precisionsKriging->computeQuadratic(_workingData);
}

double SPDE::computeLogLike(int nbsimus, int seed) const
{
  if (!_isCoeffsComputed)
  {
    _computeDriftCoeffs();
  }
  return - 0.5 * (computeLogDet(nbsimus,seed) + computeQuad()) ;
}

double SPDE::computeProfiledLogLike(int nbsimus, int seed) const
{
  _isCoeffsComputed = false; // we assume that covariance parameters have changed when using this function
  //  so driftCoeffs have to be recomputed

  return computeLogLike(nbsimus,seed);
}

void SPDE::_computeDriftCoeffs() const
{
  if (!_isCoeffsComputed)
   {
    _isCoeffsComputed = true;
    if (_requireCoeffs)
    {
      _driftCoeffs = _precisionsKriging->computeCoeffs(_data->getColumnByLocator(ELoc::Z,0,true),_driftTab);
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
