
#include "API/SPDE.hpp"
#include "Model/ANoStat.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/AException.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include <iostream>
SPDE::SPDE(Model& model,const Db& field,ANoStat* nostat,const Db* dat)
{
    init(model,field,nostat,dat);
}

void SPDE::init(Model& model,const Db& field,ANoStat* nostat,const Db* dat)
{
  double nugget=0.;
  double totalSill = 0.;
  ShiftOpCs* shiftOp;
  PrecisionOpCs precision;
  MeshETurbo* mesh;


  for(int icov = 0 ; icov < model.getCovaNumber();icov++)
  {
    const auto cova = model.getCova(icov);

    if(cova->getType()==COV_NUGGET)
    {
      nugget = cova->getSill(0,0);
    }
    else if(cova->getType()==COV_BESSEL_K)
    {
      std::cout<<"Bessel"<<std::endl;
      totalSill += cova->getSill(0,0);
      mesh = createMeshing(*cova,field,14.,0.2);
      mesh->display(0);
      shiftOp = new ShiftOpCs(mesh, &model,&field, nostat);
     // delete mesh;

     _pileShiftOp.push_back(shiftOp);
      _precistionLists.push_back(PrecisionOpCs(shiftOp,cova,POPT_MINUSHALF));

      if(dat!=nullptr)
      {

      }
    }
    else
    {
      my_throw("SPDE is only implemented for Matérn covariances (BESSEL_K)");
    }
  }
  if(nugget == 0.)
  {
    nugget = 0.01 * totalSill;
  }
  _precisionsKriging.setNugget(nugget);
}

MeshETurbo* SPDE::createMeshing(const CovAniso & cova, const Db& field,double discr,double ext)
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
SPDE::~SPDE()
{
  for(auto &e : _pileShiftOp)
  {
    delete &e;
  }
}

