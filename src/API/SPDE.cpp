#include "API/SPDE.hpp"
#include "Model/ANoStat.hpp"
#include "Covariances/CovAniso.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/AException.hpp"
#include "Basic/Law.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Db/Db.hpp"
#include <iostream>

SPDE::SPDE()
{

}

SPDE::SPDE(Model& model,const Db& field,const Db* dat,ENUM_CALCUL_MODE calc)
{
    init(model,field,dat,calc);
}

SPDE::~SPDE()
{
  _purge();
}

void SPDE::_purge()
{

}

void SPDE::init(Model& model, const Db& field, const Db* dat,ENUM_CALCUL_MODE calc)
{
  _purge();
  _model = &model;
  _calcul = calc;
  _data =  dat;
  VectorDouble varianceData;
  double totalSill = 0.;
  double nugget = 0.;
  ShiftOpCs* shiftOp;
  PrecisionOpCs* precision;
  MeshETurbo* mesh;
  ProjMatrix* proj;

  for(int icov = 0 ; icov < model.getCovaNumber();icov++)
  {
    const auto cova = model.getCova(icov);

    if (cova->getType() == COV_NUGGET)
    {
      nugget = cova->getSill(0,0);
    }
    else if (cova->getType() == COV_BESSEL_K)
    {
      totalSill += cova->getSill(0, 0);
      if(_calculSimu())
      {
        mesh = _createMeshing(*cova, field, 14., 0.2);
        _simuMeshing.push_back(mesh);
        shiftOp = new ShiftOpCs(mesh, &model, &field);
        precision = new PrecisionOpCs(shiftOp, cova, POPT_MINUSHALF);
        _pileShiftOp.push_back(shiftOp);
        _pilePrecisions.push_back(precision);
        proj = new ProjMatrix(_data,mesh);
        _pileProjMatrix.push_back(proj);
        _precisionsSimu.push_back(precision,proj);
        _workingSimu.push_back(VectorDouble(shiftOp->getSize()));
      }
      if(_calculKriging())
      {
        mesh = _createMeshing(*cova, field, 11., 0.2);
        _krigingMeshing.push_back(mesh);
        shiftOp = new ShiftOpCs(mesh, &model, &field);
        precision = new PrecisionOpCs(shiftOp, cova, POPT_ONE);
        proj = new ProjMatrix(_data,mesh);
        _pileShiftOp.push_back(shiftOp);
        _pilePrecisions.push_back(precision);
        _pileProjMatrix.push_back(proj);
        _precisionsKriging.push_back(precision,proj);
        _workKriging.push_back(VectorDouble(shiftOp->getSize()));
      }
    }
    else
    {
      my_throw("SPDE is only implemented for Matérn covariances (BESSEL_K)");
    }
  }

  // Evaluation of the variance at data point (nugget + measurement error or minimum proportion of total sill)
  if (_calculKriging())
  {
    if(dat->getVarianceErrorNumber()>0)
    {
      varianceData = dat->getFieldByLocator(LOC_V,0,true);
      for (int iech = 0; iech < dat->getActiveSampleNumber(); iech++)
      {
        double *temp = &varianceData[iech];
        *temp = MAX(*temp+nugget,0.01 * totalSill);
      }
    }
    else
    {
      varianceData = VectorDouble(dat->getActiveSampleNumber());
      for (int iech = 0; iech < dat->getActiveSampleNumber(); iech++)
      {
        varianceData[iech] = MAX(nugget, 0.01 * totalSill);
      }
    }
  }
  _precisionsKriging.setVarianceData(varianceData);
  _precisionsSimu.setVarianceData(varianceData);
}

void SPDE::computeKriging(const VectorDouble& datVect) const
{
  VectorVectorDouble rhs = _precisionsKriging.computeRhs(datVect);
  _precisionsKriging.evalInverse(rhs,_workKriging);
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
      gauss = ut_vector_simulate_gaussian(_simuMeshing[icov]->getNApices());
      _precisionsSimu.simulateOnMeshing(gauss,_workingSimu);
    }
 }
}

void SPDE::computeSimuCond(int nbsimus, int seed) const
{
  computeSimuNonCond(nbsimus,seed);
  VectorDouble temp(_data->getActiveSampleNumber());
  _precisionsSimu.simulateOnDataPointFromMeshings(_workingSimu,temp);
  ut_vector_multiply_inplace(temp,-1.);
  ut_vector_add_inplace(temp,_data->getFieldByLocator(LOC_Z,0,true));
  computeKriging(temp);
}

void SPDE::compute(int nbsimus, int seed) const
{
  if(_calcul == CALCUL_KRIGING)
  {
    computeKriging(_data->getFieldByLocator(LOC_Z,0,true));
  }

  if(_calcul == CALCUL_SIMUNONCOND)
  {
    computeSimuNonCond(nbsimus,seed);

  }
  if(_calcul == CALCUL_SIMUCOND)
  {
    computeSimuCond(nbsimus,seed);

  }
}

MeshETurbo* SPDE::_createMeshing(const CovAniso & cova,
                                const Db& field,
                                double discr,
                                double ext)
{
  VectorDouble extendMin,extendMax;
  int dim = cova.getNDim();
  for(int idim = 0;idim<(int)dim;idim++)
  {
      auto limits = field.getExtrema(idim,true);
      extendMin.push_back(limits[0]);
      extendMax.push_back(limits[1]);
  }
  VectorDouble cellSize = VectorDouble(dim,ut_vector_min(cova.getRanges())/discr);

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
  mesh->display();
  return mesh;
}

void SPDE::query(Db* db, NamingConvention namconv)
{
  VectorDouble temp(db->getActiveSampleNumber());
  VectorDouble result(db->getActiveSampleNumber(),0.);
  String suffix;
  if(_calcul == CALCUL_KRIGING)
  {
    for(int i = 0 ; i< (int)_krigingMeshing.size(); i++)
    {
      ProjMatrix proj(db,_krigingMeshing[i]);
      proj.mesh2point(_workKriging[i],temp);
      ut_vector_add_inplace(result,temp);
    }
    suffix = "kriging";
  }
  else if(_calcul == CALCUL_SIMUNONCOND)
  {
    for(int i = 0 ; i< (int)_simuMeshing.size(); i++)
    {
      ProjMatrix proj(db,_simuMeshing[i]);
      proj.mesh2point(_workingSimu[i],temp);
      // TODO add nugget
      ut_vector_add_inplace(result,temp);
    }
    suffix = "simu";
  }
  else if(_calcul == CALCUL_SIMUCOND)
  {
    for(int i = 0 ; i< (int)_simuMeshing.size(); i++)
    {
      ProjMatrix projSimu(db,_simuMeshing[i]);
      projSimu.mesh2point(_workingSimu[i],temp);
      ut_vector_add_inplace(result,temp);
      ProjMatrix projKriging(db,_krigingMeshing[i]);
      projKriging.mesh2point(_workKriging[i],temp);
      ut_vector_add_inplace(result,temp);
      // TODO add nugget
     }
      suffix = "condSimu";
  }

  int iptr = db->addFields(result,"SPDE",LOC_Z,true,TEST);
  namconv.setNamesAndLocators(_data,LOC_Z,1,db,iptr,suffix,1,true);
}

VectorDouble SPDE::computeCoeffs(const VectorVectorDouble& x)const
{
  return _precisionsKriging.computeCoeffs(_data->getFieldByLocator(LOC_Z,0,true),x);
}
