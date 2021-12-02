/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/*                                                                            */
/* This file is meant to demonstrate the process of performing                */
/* non-conditional simulations using in turn simpgs and simbipgs              */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"
#include "Db/Db.hpp"
#include "Covariances/ECov.hpp"
#include "Covariances/CovAniso.hpp"
#include "Variogram/Vario.hpp"
#include "Neigh/Neigh.hpp"
#include "Model/Model.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleShift.hpp"
#include "LithoRule/RuleShadow.hpp"
#include "LithoRule/RuleProp.hpp"


/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])

{
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("simbiPGS-");
  int error = 0;
  int ndim = 2;
  int nbsimu = 2;
  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  CovContext ctxt(1,2,1.); // use default space

  // Prepare the Discrete process with Discretized Option
  set_test_discrete(false);

  // Creating an output Grid Db
  Db dbgrid= Db({100,100},{0.01,0.01},{0.,0.});

  // Creating the proportions for simPGS
  VectorDouble props1({0.2, 0.5, 0.3});

  // Creating the Model(s) of the Underlying GRF(s)
  Model model1(ctxt);
  double range1 = 0.2;
  CovAniso cova1(ECov::BESSEL_K,range1,1.,1.,ctxt);
  model1.addCova(&cova1);
  model1.display();
  model1.serialize("PGSmodel1.ascii");

  Model model2(ctxt);
  double range2 = 0.3;
  CovAniso cova2(ECov::EXPONENTIAL,range2,1.,1.,ctxt);
  model2.addCova(&cova2);
  model2.display();
  model2.serialize("PGSmodel2.ascii");

  Model model3(ctxt);
  double range3 = 0.2;
  CovAniso cova3(ECov::BESSEL_K,range3,1.,1.,ctxt);
  model3.addCova(&cova3);
  model3.display();
  model3.serialize("PGSmodel3.ascii");

  Model model4(ctxt);
  double range4 = 0.1;
  CovAniso cova4(ECov::SPHERICAL,range4,1.,1.,ctxt);
  model4.addCova(&cova4);
  model4.display();
  model4.serialize("PGSmodel4.ascii");

  // Creating the Neighborhood
  Neigh neigh = Neigh(ndim);
  neigh.display();

  // Creating the Rules
  Rule rule1({"S","S","F1","F2","F3"});
  rule1.display();
  rule1.serialize("PGSrule1.ascii");

  // Creating the RuleProp structure for simPGS
  RuleProp ruleprop1 = RuleProp(&rule1, props1);

  // Perform a non-conditional PGS simulation on a grid
  error = simpgs(nullptr,&dbgrid,&ruleprop1,&model1,&model2,&neigh,nbsimu);
  dbgrid.setName(ELoc::FACIES,"PGS-Facies");
  dbgrid.display();
  dbgrid.serialize("simupgs.ascii");

  // Creating the RuleProp for simBiPGS
  VectorDouble props2({0.1, 0.2, 0.1, 0.3, 0.1, 0.2});
  Rule rule2({"S","F1","F2"});
  rule2.display();
  rule2.serialize("PGSrule2.ascii");
  RuleProp rulepropbi = RuleProp(&rule1, &rule2, props2);

  // Perform a non-conditional BiPGS simulation on a grid
  error = simbipgs(nullptr,&dbgrid,&rulepropbi,
                   &model1,&model2,&model3,&model4,&neigh,nbsimu);
  dbgrid.setName(ELoc::FACIES,"BiPGS-Facies");
  dbgrid.display();
  dbgrid.serialize("simubipgs.ascii");

  // Performing a PGS simulation using Shift
  VectorDouble shift = {0.2, 0.3};
  VectorDouble propshift = { 0.1, 0.2, 0.3, 0.4 };
  RuleShift ruleshift({"S","S","S","F1","F2","F3","F4"},shift);
  ruleshift.display(1);
  ruleshift.serialize("PGSruleshift.ascii");

  RuleProp rulepropshift = RuleProp(&ruleshift, propshift);

  // Perform a non-conditional PGS Shift simulation on a grid
  error = simpgs(nullptr,&dbgrid,&rulepropshift,&model1,nullptr,&neigh,nbsimu);
  dbgrid.setName(ELoc::FACIES,"PGS-Shift-Facies");
  dbgrid.display();
  dbgrid.serialize("simushiftpgs.ascii");

  // Performing a PGS simulation using Shadow
  double slope = 0.5;
  double shdown = -0.2;
  double shdsup = +0.5;
  RuleShadow ruleshadow = RuleShadow(slope,shdsup,shdown,shift);
  ruleshadow.display(1);
  ruleshadow.serialize("PGSruleshadow.ascii");

  VectorDouble propshadow = { 0.4, 0.2, 0.3 };
  RuleProp rulepropshadow = RuleProp(&ruleshadow, propshadow);

  // Perform a non-conditional PGS Shadow simulation on a grid
  error = simpgs(nullptr,&dbgrid,&rulepropshadow,&model1,nullptr,&neigh,nbsimu);
  dbgrid.setName(ELoc::FACIES,"PGS-Shadow-Facies");
  dbgrid.display();
  dbgrid.serialize("simushadowpgs.ascii");

  return(error);
}
