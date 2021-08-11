
#include "API/SPDE.hpp"
#include "Model/ANoStat.hpp"
#include "Covariances/CovAniso.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/AException.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "Db/Db.hpp"
#include <iostream>

SPDE::SPDE(Model& model,const Db& field,const Db* dat)
{
    init(model,field,dat);
}

SPDE::~SPDE()
{
  for(auto &e : _pileShiftOp)
  {
    //delete &e;
  }
}

void SPDE::init(Model& model, const Db& field, const Db* dat)
{
  _model = &model;
  VectorDouble varianceData;
  double totalSill = 0.;
  double nugget = 0.;
  ShiftOpCs* shiftOp;
  PrecisionOpCs precision;
  MeshETurbo* mesh;

  for(int icov = 0 ; icov < model.getCovaNumber();icov++)
  {
    const auto cova = model.getCova(icov);

    if (cova->getType() == COV_NUGGET)
    {
      nugget = cova->getSill(0,0);
    }
    else if (cova->getType() == COV_BESSEL_K)
    {
      std::cout << "Bessel" << std::endl;
      totalSill += cova->getSill(0, 0);
      mesh = createMeshing(*cova, field, 14., 0.2);
      mesh->display(0);
      shiftOp = new ShiftOpCs(mesh, &model, &field);
      // delete mesh;

      _pileShiftOp.push_back(shiftOp);
      _precistionLists.push_back(PrecisionOpCs(shiftOp, cova, POPT_MINUSHALF));
    }
    else
    {
      my_throw("SPDE is only implemented for Matérn covariances (BESSEL_K)");
    }
  }

  if (dat != nullptr)
  {
    if(dat->getVarianceErrorNumber()>0)
    {
      varianceData = dat->getFieldByLocator(LOC_V,0,true);
    }
    for (int iech = 0; iech < dat->getActiveSampleNumber(); iech++)
    {
      double *temp = &varianceData[iech];
      *temp = MAX(*temp+nugget,0.01 * totalSill);
    }
  }
  _precisionsKriging.setVarianceData(varianceData);
}

MeshETurbo* SPDE::createMeshing(const CovAniso & cova,
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
