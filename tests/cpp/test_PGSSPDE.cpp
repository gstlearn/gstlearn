/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/*                                                                            */
/* This file is meant to demonstrate the process of using PGS                 */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"

#include "Enum/ECov.hpp"

#include "Model/Model.hpp"
#include "API/SPDE.hpp"
#include "API/PGSSPDE.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Basic/String.hpp"
#include "Basic/File.hpp"
#include "Basic/VectorHelper.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleProp.hpp"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix acovalistlgebra
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("PGSSPDE-");
  int error = 0;
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  CovContext ctxt(1,2,1.); // use default space

  // Prepare the Discrete process with Discretized Option
  set_test_discrete(false);

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 10;
  Db* db = Db::createFromBox(nech,{0.,0.},{1.,1.}, 43431);
  db->display(); // TODO : please use FLAG_STATS only when available

  auto nx={ 101,101 };
  DbGrid* workingDbc = DbGrid::create(nx);
  DbGrid* dbprop = DbGrid::create({100,100},{0.01,0.01});

  VectorDouble props({0.2, 0.5, 0.3});
  int nfac = props.size();
  VectorString names = generateMultipleNames("Props",nfac);
  for (int ifac = 0; ifac < nfac; ifac++)
    dbprop->addColumnsByConstant(1,props[ifac],names[ifac]);
  dbprop->setLocators(names,ELoc::P);

  // Creating the Model(s) of the Underlying GRF(s)
  Model model1(ctxt);
  CovLMC covs1(ctxt.getSpace());
  double range1 = 20;
  CovAniso cova1(ECov::BESSEL_K,range1,1.,1.,ctxt);
  covs1.addCov(&cova1);
  model1.setCovList(&covs1);
  model1.display();
  (void) model1.dumpToNF("truemodel1.ascii");

  Model model2(ctxt);
  CovLMC covs2(ctxt.getSpace());
  double range2 = 40;
  CovAniso cova2(ECov::BESSEL_K,range2,2.,1.,ctxt);
  covs2.addCov(&cova2);
  model2.setCovList(&covs2);
  model2.display();
  (void) model2.dumpToNF("truemodel2.ascii");

  std::vector<Model*> models;
  models.push_back(&model1);
  models.push_back(&model2);

  // Creating the Rule
  Rule* rule = Rule::createFromNames({"S","T","F1","F2","F3"});
  RuleProp* ruleprop = RuleProp::createFromRule(rule, props);

  auto ndata = 100;
  Db* dat = Db::createFromBox(ndata, {0.,0.}, {100.,100.}, 32432);
  VectorDouble z = VH::simulateGaussian(ndata);
  dat->addColumns(z,"variable",ELoc::Z);

  PGSSPDE sCond(models,workingDbc,ruleprop,dat);
  PGSSPDE sNonCond(models,workingDbc,ruleprop);

  PGSSPDE* spgs = &sNonCond;
  spgs->simulate();
  spgs->query(workingDbc);
  DbStringFormat dbfmt(FLAG_STATS,{"Facies"});
  workingDbc->display(&dbfmt);
  (void) workingDbc->dumpToNF("pgs.ascii");

  delete db;
  delete workingDbc;
  delete dbprop;
  delete dat;
  delete rule;
  delete ruleprop;
  return(error);
}
