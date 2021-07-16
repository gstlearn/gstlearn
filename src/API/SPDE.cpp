
#include "API/SPDE.hpp"
#include "Model/ANoStat.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/AException.hpp"
#include "LinearOp/ShiftOpCs.hpp"

SPDE::SPDE(Model& model,const Db& field,ANoStat* nostat,const Db* dat)
{
  init(model,field,nostat,dat);

}

void SPDE::init(Model& model,const Db& field,ANoStat* nostat,const Db* dat)
{
  double nugget=0.;
  double totalSill = 0.;
  MeshETurbo mesh;
  ShiftOpCs shiftOp;

  for(int icov = 0 ; icov < model.getCovaNumber();icov++)
  {
    auto cova = model.getCova(icov);
    if(cova->getType()==COV_NUGGET)
    {
      nugget = cova->getSill(0,0);
    }
    else if(cova->getType()==COV_BESSEL_K)
    {
      totalSill += cova->getSill(0,0);

      mesh = createMeshing(*cova,field,14);
      shiftOp.initFromMesh(&mesh,&model,&field,nostat,0,icov);

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

MeshETurbo SPDE::createMeshing(const CovAniso & cova, const Db& field,double discr)
{
  VectorDouble extendMin,extendMax;
  int dim = cova.getNDim();
  for(int idim = 0;idim<(int)dim;idim++)
  {
      auto limits = field.getExtrema(idim,true);
      extendMin.push_back(limits[0]);
      extendMax.push_back(limits[1]);
  }
  VectorDouble cellSize = VectorDouble(dim,ut_vector_min(cova.getAnisoCoeffs())/discr);

  MeshETurbo mesh;
  mesh.initFromExtend(extendMin,extendMax,cellSize,field.getRotMat());
  return mesh;
}
SPDE::~SPDE()
{
  // TODO Auto-generated destructor stub
}

