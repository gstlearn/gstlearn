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
<<<<<<< HEAD
  PGSSPDE(std::vector<Model*> models,
         const Db& field,
         RuleProp ruleprop,
         const Db* dat=nullptr);
  void simulate(int seed= 32145,int nitergibbs = 0) const;
  void simulateNonCond(int seed = 32145) const;
=======
  PGSSPDE(std::vector<Model> models,
          const Db& field,
          RuleProp ruleprop,
          const Db* dat = nullptr);
>>>>>>> 219703517e4414ed713ae3003b86fd336c9c5573
  void gibbs(int niter) const;
  void query(Db* db,bool keepGauss=false) const;
  virtual ~PGSSPDE();
private:
  std::vector<SPDE> _spdeTab;
  RuleProp _ruleProp;
  mutable Db* _workingDb;
  ENUM_CALCUL_MODE _calcul;
};
