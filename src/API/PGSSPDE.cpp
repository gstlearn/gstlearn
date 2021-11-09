#include "API/PGSSPDE.hpp"

PGSSPDE::PGSSPDE(std::vector<Model*> models,
                 const Db& field,
                 RuleProp ruleprop,
                 const Db* dat)
: _data()
, _spdeTab()
, _ruleProp()
, _calcul()
{
  _calcul = (dat == nullptr) ? ESPDECalcMode::SIMUNONCOND : ESPDECalcMode::SIMUCOND;
  for(auto &e : models)
  {
    _spdeTab.push_back(SPDE(*e,field,dat,_calcul));
  }
  _ruleProp = ruleprop;
}

void PGSSPDE::simulate(int seed,int /*nitergibbs*/) const
{
  if(_calcul==ESPDECalcMode::SIMUNONCOND)
  {
    simulateNonCond(seed);
  }
  if(_calcul==ESPDECalcMode::SIMUCOND)
  {
    gibbs(1);
  }
}

void PGSSPDE::simulateNonCond(int seed) const
{
  for(int i = 0; i < (int)_spdeTab.size();i++)
  {
    int curseed = i==0? seed:0;
    _spdeTab[i].compute(1,curseed);
  }
}

void PGSSPDE::query(Db* db,bool keepGauss) const
{

  int ngrf = (int)_spdeTab.size();
  VectorString names = generateMultipleNames("simuGauss",ngrf);

  for(int i = 0; i < ngrf;i++)
  {
   int iptr = _spdeTab[i].query(db,NamingConvention("simGauss"));
   db->setName(iptr,names[i]);
  }

  db->setLocator(names,ELoc::Z);
  db->display();
  _ruleProp.gaussToCategory(db,NamingConvention("categories"));

  if(!keepGauss)
  {
    db->deleteField(names);
  }
}

void PGSSPDE::gibbs(int /*niter*/) const
{
   _ruleProp.categoryToThresh(_data);
}

PGSSPDE::~PGSSPDE()
{
  // TODO Auto-generated destructor stub
}

