#include "API/SPDE.hpp"
#include "Model/ANoStat.hpp"
#include "Covariances/CovAniso.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/AException.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Db/Db.hpp"
#include <iostream>

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
        shiftOp = new ShiftOpCs(mesh, &model, &field);
        precision = new PrecisionOpCs(shiftOp, cova, POPT_MINUSHALF);
        _pileShiftOp.push_back(shiftOp);
        _pilePrecisions.push_back(precision);
      }
      if(_calculKriging())
      {
        mesh = _createMeshing(*cova, field, 11., 0.2);
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

  if (dat != nullptr && _calculKriging())
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
}

void SPDE::computeKriging() const
{
  VectorDouble datVect = _data->getFieldByLocator(LOC_Z,0,true);
  VectorVectorDouble rhs = _precisionsKriging.computeRhs(datVect);
  _precisionsKriging.evalInverse(rhs,_workKriging);
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
  return new MeshETurbo(nx,cellSize,x0,field.getRotMat());
}
