#include <API/PGSSPDE.hpp>
#include "geoslib_enum.h"

PGSSPDE::PGSSPDE(std::vector<Model> models,
         const Db& field,
         RuleProp ruleprop,
         const Db* dat)
{
  ENUM_CALCUL_MODE calc = (dat == nullptr) ? CALCUL_SIMUNONCOND:CALCUL_SIMUCOND;
  for(auto &e : models)
  {
    _spdeTab.push_back(SPDE(e,field,dat,calc));
  }

  _ruleProp = ruleprop;
}


void PGSSPDE::gibbs(int niter) const
{

}

PGSSPDE::~PGSSPDE()
{
  // TODO Auto-generated destructor stub
}

