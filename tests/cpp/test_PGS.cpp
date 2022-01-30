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
/* This file is meant to demonstrate the process of using PGS                 */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"
#include "Variogram/Vario.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Model/Model.hpp"
#include "LithoRule/RuleProp.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/String.hpp"
#include "Covariances/ECov.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])

{
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("PGS-");
  int error = 0;
  int ndim  = 2;
  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  CovContext ctxt(1,2,1.); // use default space

  // Prepare the Discrete process with Discretized Option
  set_test_discrete(false);
  bool flagStationary = false;

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 1000;
  Db* db = Db::createFromBox(nech,{0.,0.},{1.,1.});
  DbStringFormat dbfmt;
  dbfmt.setParams(FLAG_STATS);
  db->display(&dbfmt);

  Dbgrid* dbprop = Dbgrid::create({100,100},{0.01,0.01});

  VectorDouble props({0.2, 0.5, 0.3});
  int nfac = props.size();
  VectorString names = generateMultipleNames("Props",nfac);
  for (int ifac = 0; ifac < nfac; ifac++)
    dbprop->addFieldsByConstant(1,props[ifac],names[ifac],ELoc::P,ifac);
  dbprop->display();

  // Creating the Model(s) of the Underlying GRF(s)
  Model model1(ctxt);
  CovLMC covs1(ctxt.getSpace());
  double range1 = 0.2;
  CovAniso cova1(ECov::BESSEL_K,range1,1.,1.,ctxt);
  covs1.addCov(&cova1);
  model1.setCovList(&covs1);
  model1.display();
  (void) model1.dumpToNF("truemodel1.ascii");

  Model model2(ctxt);
  CovLMC covs2(ctxt.getSpace());
  double range2 = 0.3;
  CovAniso cova2(ECov::EXPONENTIAL,range2,1.,1.,ctxt);
  covs2.addCov(&cova2);
  model2.setCovList(&covs2);
  model2.display();
  (void) model2.dumpToNF("truemodel2.ascii");

  // Creating the Neighborhood
  NeighUnique* neighU = NeighUnique::create(ndim, false);
  neighU->display();

  // Creating the Rule
  Rule* rule = Rule::createFromNames({"S","T","F1","F2","F3"});
  rule->display();
  rule->dumpToNF("truerule.ascii");
  RuleProp* ruleprop;
  if (flagStationary)
    ruleprop = RuleProp::createFromRule(rule, props);
  else
    ruleprop = RuleProp::createFromRuleAndDb(rule, dbprop);

  // Perform a non-conditional simulation on the Db
  error = simpgs(nullptr,db,ruleprop,&model1,&model2,neighU);
  db->setLocator(db->getLastName(),ELoc::Z);
  (void) db->dumpToNF("simupgs.ascii");

  // Design of several VarioParams
  int nlag1 = 19;
  DirParam dirparam1 = DirParam(ndim, nlag1, 0.5 / nlag1);
  VarioParam varioparam1;
  varioparam1.addDirs(dirparam1);

  int nlag2 = 3;
  DirParam dirparam2 = DirParam(ndim, nlag2, 0.1 );
  VarioParam varioparam2;
  varioparam2.addDirs(dirparam2);

  // Determination of the variogram of the Underlying GRF
  Vario* vario = variogram_pgs(db,&varioparam1,ruleprop);
  vario->display();
  Vario vario1 = Vario(*vario);
  vario1.reduce({0},VectorInt(),true);
  Vario vario2 = Vario(*vario);
  vario2.reduce({1},VectorInt(),true);
  vario1.display();
  vario2.display();

  // Fitting the experimental variogram of Underlying GRF (with constraint that total sill is 1)
  Model modelPGS1(ctxt);
  Model modelPGS2(ctxt);
  Option_AutoFit option = Option_AutoFit();
  option.setConstantSillValue(1.);

  std::vector<ECov> covs {ECov::BESSEL_K, ECov::EXPONENTIAL};
  modelPGS1.fit(&vario1,covs,true,option);
  modelPGS1.display();

  (void) vario1.dumpToNF("variopgs1.ascii");
  (void) modelPGS1.dumpToNF("modelfitpgs1.ascii");

  modelPGS2.fit(&vario2,covs,true,option);
  modelPGS2.display();

  (void) vario2.dumpToNF("variopgs2.ascii");
  (void) modelPGS2.dumpToNF("modelfitpgs2.ascii");

  RuleProp* ruleprop2;
  if (flagStationary)
    ruleprop2 = RuleProp::createFromRule((Rule*) NULL, props);
  else
    ruleprop2 = RuleProp::createFromDb(dbprop, VectorDouble());
  error = ruleprop2->fit(db, &varioparam2, 2, true);
  ruleprop2->getRule()->display();
  (void) ruleprop2->getRule()->dumpToNF("ruleFit.ascii");

  modelPGS1.display();
  Vario* varioDerived = model_pgs(db, &varioparam1, ruleprop2, &modelPGS1, &modelPGS2);
  modelPGS1.display();
  varioDerived->dumpToNF("modelpgs.ascii");
  varioDerived->display();

  Vario varioIndic = Vario(&varioparam1, db);
  varioIndic.computeIndic("vg");
  (void) varioIndic.dumpToNF("varioindic.ascii");

  modelPGS1.display();

  delete db;
  delete dbprop;
  delete neighU;
  delete rule;
  delete ruleprop;
  delete ruleprop2;
  return(error);
}
