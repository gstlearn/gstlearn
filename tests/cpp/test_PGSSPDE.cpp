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
/* This file is meant to demonstrate the process of using PGS                 */
/*                                                                            */
/******************************************************************************/
#include "Basic/Law.hpp"
#include "geoslib_f.h"

#include "Enum/ECov.hpp"

#include "Model/Model.hpp"
#include "API/SPDE.hpp"
#include "API/PGSSPDE.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
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
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("test_PGSSPDE-");
  int error = 0;
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Prepare the Discrete process with Discretized Option
  set_test_discrete(false);

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 10;
  Db* db = Db::createFromBox(nech,{0.,0.},{1.,1.}, 43431);
  db->display(); // TODO : please use FLAG_STATS only when available

  auto nx={ 101,101 };
  DbGrid* grid = DbGrid::create(nx);
  DbGrid* dbprop = DbGrid::create({100,100},{0.01,0.01});

  VectorDouble props({0.2, 0.5, 0.3});
  int nfac = (int) props.size();
  VectorString names = generateMultipleNames("Props",nfac);
  for (int ifac = 0; ifac < nfac; ifac++)
    dbprop->addColumnsByConstant(1,props[ifac],names[ifac]);
  dbprop->setLocators(names,ELoc::P);

  // Creating the Model(s) of the Underlying GRF(s)
  double range1 = 20;
  Model* model1 = Model::createFromParam(ECov::MATERN, range1, 1., 1.);
  model1->display();

  double range2 = 40;
  Model* model2 = Model::createFromParam(ECov::MATERN, range2, 1., 2.);
  model2->display();

  std::vector<Model*> models;
  models.push_back(model1);
  models.push_back(model2);

  // Creating the Rule
  Rule* rule = Rule::createFromNames({"S","T","F1","F2","F3"});
  RuleProp* ruleprop = RuleProp::createFromRule(rule, props);
  rule->dumpToNF("rule.ascii");

  auto ndata = 100;
  Db* dat = Db::createFromBox(ndata, {0.,0.}, {100.,100.}, 32432);
  VectorDouble z = VH::simulateGaussian(ndata);
  dat->addColumns(z,"variable",ELoc::Z);

  PGSSPDE sNonCond(models,grid,ruleprop);
  law_set_random_seed(133672);
  sNonCond.compute(grid, 0, NamingConvention("Facies-NC"));

  PGSSPDE sCond(models,grid,ruleprop,dat);
  law_set_random_seed(133272);
  sCond.compute(grid, 0, NamingConvention("Facies-CD"));

  DbStringFormat dbfmt(FLAG_STATS,{"Facies"});
  grid->display(&dbfmt);
  (void) grid->dumpToNF("pgs.ascii");

  delete db;
  delete grid;
  delete dbprop;
  delete dat;
  delete rule;
  delete ruleprop;
  return(error);
}
