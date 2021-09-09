#pragma once

#include "Basic/NamingConvention.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "API/SPDE.hpp"
#include "LithoRule/RuleProp.hpp"
#include "geoslib_enum.h"
#include <vector>

class PGSSPDE
{
public:
  PGSSPDE(std::vector<Model*> models,
         const Db& field,
         RuleProp ruleprop,
         const Db* dat=nullptr);
  void simulate(int seed= 32145,int nitergibbs = 0) const;
  void simulateNonCond(int seed = 32145) const;
  void gibbs(int niter) const;
  void query(Db* db,bool keepGauss=false) const;
  virtual ~PGSSPDE();
private:
  Db* _data;
  std::vector<SPDE> _spdeTab;
  RuleProp _ruleProp;
  ENUM_CALCUL_MODE _calcul;
};
