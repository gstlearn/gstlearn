/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/*                                                                            */
/* This file is meant to demonstrate the process of performing                */
/* non-conditional simulations using in turn simpgs and simbipgs              */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"

#include "Enum/ECov.hpp"

#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Variogram/Vario.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleShift.hpp"
#include "LithoRule/RuleShadow.hpp"
#include "LithoRule/RuleProp.hpp"

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("simbiPGS-");
  int error = 0;
  int ndim = 2;
  int nbsimu = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  CovContext ctxt(1,2,1.); // use default space
  DbStringFormat dbfmt;

  // Prepare the Discrete process with Discretized Option
  set_test_discrete(false);

  // Creating an output Grid Db
  DbGrid* dbgrid = DbGrid::create({100,100},{0.01,0.01},{0.,0.});

  // Creating the proportions for simPGS
  VectorDouble props1({0.2, 0.5, 0.3});

  // Creating the Model(s) of the Underlying GRF(s)
  Model model1(ctxt);
  CovLMC covs1(ctxt.getSpace());
  double range1 = 0.2;
  CovAniso cova1(ECov::BESSEL_K,range1,1.,1.,ctxt);
  covs1.addCov(&cova1);
  model1.setCovList(&covs1);
  model1.display();
  (void) model1.dumpToNF("PGSmodel1.ascii");

  Model model2(ctxt);
  CovLMC covs2(ctxt.getSpace());
  double range2 = 0.3;
  CovAniso cova2(ECov::EXPONENTIAL,range2,1.,1.,ctxt);
  covs2.addCov(&cova2);
  model2.setCovList(&covs2);
  model2.display();
  (void) model2.dumpToNF("PGSmodel2.ascii");

  Model model3(ctxt);
  CovLMC covs3(ctxt.getSpace());
  double range3 = 0.2;
  CovAniso cova3(ECov::BESSEL_K,range3,1.,1.,ctxt);
  covs3.addCov(&cova3);
  model3.setCovList(&covs3);
  model3.display();
  (void) model3.dumpToNF("PGSmodel3.ascii");

  Model model4(ctxt);
  CovLMC covs4(ctxt.getSpace());
  double range4 = 0.1;
  CovAniso cova4(ECov::SPHERICAL,range4,1.,1.,ctxt);
  covs4.addCov(&cova4);
  model4.setCovList(&covs4);
  model4.display();
  (void) model4.dumpToNF("PGSmodel4.ascii");

  // Creating the Neighborhood
  NeighUnique* neighU = NeighUnique::create();
  neighU->display();

  // Creating the Rules
  Rule* rule1 = Rule::createFromNames({"S","S","F1","F2","F3"});
  rule1->display();
  (void) rule1->dumpToNF("PGSrule1.ascii");

  // Creating the RuleProp structure for simPGS
  RuleProp* ruleprop1 = RuleProp::createFromRule(rule1, props1);

  // Perform a non-conditional PGS simulation on a grid
  error = simpgs(nullptr,dbgrid,ruleprop1,&model1,&model2,neighU,nbsimu);
  dbgrid->setNameByLocator(ELoc::FACIES,"PGS-Facies");
  dbfmt = DbStringFormat(FLAG_STATS,{"PGS-Facies*"});
  dbgrid->display(&dbfmt);
  (void) dbgrid->dumpToNF("simupgs.ascii");

  // Creating the RuleProp for simBiPGS
  VectorDouble props2({0.1, 0.2, 0.1, 0.3, 0.1, 0.2});
  Rule* rule2 = Rule::createFromNames({"S","F1","F2"});
  rule2->display();
  (void) rule2->dumpToNF("PGSrule2.ascii");
  RuleProp* rulepropbi = RuleProp::createFromRules(rule1, rule2, props2);

  // Perform a non-conditional BiPGS simulation on a grid
  error = simbipgs(nullptr,dbgrid,rulepropbi,
                   &model1,&model2,&model3,&model4,neighU,nbsimu);
  dbgrid->setNameByLocator(ELoc::FACIES,"BiPGS-Facies");
  dbfmt = DbStringFormat(FLAG_STATS,{"BiPGS-Facies*"});
  dbgrid->display(&dbfmt);
  (void) dbgrid->dumpToNF("simubipgs.ascii");

  // Performing a PGS simulation using Shift
  VectorDouble shift = {0.2, 0.3};
  VectorDouble propshift = { 0.1, 0.2, 0.3, 0.4 };
  RuleShift* ruleshift = RuleShift::createFromNames({"S","S","S","F1","F2","F3","F4"},shift);
  ruleshift->display();
  (void) ruleshift->dumpToNF("PGSruleshift.ascii");

  RuleProp* rulepropshift = RuleProp::createFromRule(ruleshift, propshift);

  // Perform a non-conditional PGS Shift simulation on a grid
  error = simpgs(nullptr,dbgrid,rulepropshift,&model1,nullptr,neighU,nbsimu);
  dbgrid->setNameByLocator(ELoc::FACIES,"PGS-Shift-Facies");
  dbfmt = DbStringFormat(FLAG_STATS,{"PGS-Shift-Facies*"});
  dbgrid->display(&dbfmt);
  (void) dbgrid->dumpToNF("simushiftpgs.ascii");

  // Performing a PGS simulation using Shadow
  double slope = 0.5;
  double shdown = -0.2;
  double shdsup = +0.5;
  RuleShadow* ruleshadow = new RuleShadow(slope,shdsup,shdown,shift);
  ruleshadow->display();
  (void) ruleshadow->dumpToNF("PGSruleshadow.ascii");

  VectorDouble propshadow = { 0.4, 0.2, 0.3 };
  RuleProp* rulepropshadow = RuleProp::createFromRule(ruleshadow, propshadow);

  // Perform a non-conditional PGS Shadow simulation on a grid
  error = simpgs(nullptr,dbgrid,rulepropshadow,&model1,nullptr,neighU,nbsimu);
  dbgrid->setNameByLocator(ELoc::FACIES,"PGS-Shadow-Facies");
  dbfmt = DbStringFormat(FLAG_STATS,{"PGS-Shadow-Facies*"});
  dbgrid->display(&dbfmt);
  (void) dbgrid->dumpToNF("simushadowpgs.ascii");

  delete dbgrid;
  delete rule1;
  delete rule2;
  delete ruleshift;
  delete ruleshadow;
  delete ruleprop1;
  delete rulepropbi;
  delete rulepropshift;
  delete rulepropshadow;
  delete neighU;
  return(error);
}
