#pragma once

#include "Basic/NamingConvention.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "API/SPDE.hpp"
#include "LithoRule/RuleProp.hpp"
#include <vector>

class PGSSPDE
{
public:
  PGSSPDE(std::vector<Model> models,
         const Db& field,
         RuleProp ruleprop,
         const Db* dat=nullptr);
  void gibbs(int niter) const;
  virtual ~PGSSPDE();
private:
  std::vector<SPDE> _spdeTab;
  RuleProp _ruleProp;
};

