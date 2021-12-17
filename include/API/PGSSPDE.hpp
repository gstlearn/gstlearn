#pragma once

#include "gstlearn_export.hpp"

#include "API/ESPDECalcMode.hpp"
#include "Covariances/ECalcMember.hpp"

#include "API/SPDE.hpp"
#include "LithoRule/RuleProp.hpp"

#include <vector>

class Db;
class Model;

class GSTLEARN_EXPORT PGSSPDE
{
public:
  PGSSPDE(std::vector<Model*> models,
          const Db& field,
          RuleProp ruleprop,
          const Db* dat=nullptr);
  PGSSPDE(const PGSSPDE& r) = delete;
  PGSSPDE& operator=(const PGSSPDE& r) = delete;
  virtual ~PGSSPDE();
  void simulate(int seed= 32145,int nitergibbs = 0) const;
  void simulateNonCond(int seed = 32145) const;
  void gibbs(int niter) const;
  void query(Db* db,bool keepGauss=false) const;

private:
  Db*               _data;
  std::vector<SPDE*> _spdeTab;
  RuleProp          _ruleProp;
  ESPDECalcMode     _calcul;
};
