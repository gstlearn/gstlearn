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
#include "geoslib_f.h"

#include "Enum/ECov.hpp"

#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Covariances/KernelMatern.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "LithoRule/RuleProp.hpp"
#include "Model/Constraints.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Variogram/Vario.hpp"

using namespace gstlrn;
/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int argc, char* argv[])
{
  bessel_set_old_style(true);
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("test_PGS-");
  Id error = 0;
  Id ndim  = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  CovContext ctxt(1, 2, 1.);

  // Prepare the Discrete process with Discretized Option
  set_test_discrete(false);
  bool flagStationary = false;

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  Id nech = 1000;
  Db* db  = Db::createFromBox(nech, {0., 0.}, {1., 1.}, 432432);
  DbStringFormat dbfmt(FLAG_STATS);
  db->display(&dbfmt);

  // Create a Grid covering the area
  DbGrid* dbprop = DbGrid::create({100, 100}, {0.01, 0.01});

  // Setting constant global proportions
  VectorDouble props({0.2, 0.5, 0.3});
  Id nfac            = static_cast<Id>(props.size());
  VectorString names = generateMultipleNames("Props", nfac);
  for (Id ifac = 0; ifac < nfac; ifac++)
    dbprop->addColumnsByConstant(1, props[ifac], names[ifac], ELoc::P, ifac);
  dbprop->display();

  // Creating the Model(s) of the Underlying GRF(s)
  double range1 = 0.2;
  Model* model1 = Model::createFromParam(ECov::MATERN, range1, 1., 1.);
  model1->display();
  (void)model1->dumpToNF("truemodel1.NF");

  double range2 = 0.3;
  Model* model2 = Model::createFromParam(ECov::EXPONENTIAL, range2, 1., 1.);
  model2->display();
  (void)model2->dumpToNF("truemodel2.NF");

  // Creating the Neighborhood
  NeighUnique* neighU = NeighUnique::create();
  neighU->display();

  // Creating the Rule
  Rule* rule = Rule::createFromNames({"S", "T", "F1", "F2", "F3"});
  rule->display();
  (void)rule->dumpToNF("truerule.NF");
  RuleProp* ruleprop;
  if (flagStationary)
    ruleprop = RuleProp::createFromRule(rule, props);
  else
    ruleprop = RuleProp::createFromRuleAndDb(rule, dbprop);

  // Perform a non-conditional simulation on the Db
  error = simpgs(nullptr, db, ruleprop, model1, model2, neighU);
  db->setLocator(db->getLastName(), ELoc::Z, 0);
  (void)db->dumpToNF("simupgs.NF");
  db->display(&dbfmt);

  // Design of several VarioParams
  Id nlag1 = 19;
  DirParam dirparam1(nlag1, 0.5 / nlag1);
  VarioParam varioparam1;
  varioparam1.addDir(dirparam1);

  Id nlag2 = 3;
  DirParam dirparam2(nlag2, 0.1);
  VarioParam varioparam2;
  varioparam2.addDir(dirparam2);

  // Determination of the variogram of the Underlying GRF
  Vario* vario = variogram_pgs(db, &varioparam1, ruleprop);
  vario->display();
  Vario vario1(*vario);
  vario1.resetReduce({0}, VectorInt(), true);
  Vario vario2(*vario);
  vario2.resetReduce({1}, VectorInt(), true);
  vario1.display();
  vario2.display();

  // Fitting the experimental variogram of Underlying GRF (with constraint that total sill is 1)
  Model modelPGS1(ctxt);
  Model modelPGS2(ctxt);
  Constraints constraints;
  constraints.setConstantSillValue(1.);

  VectorECov covs {ECov::MATERN, ECov::EXPONENTIAL};
  modelPGS1.fit(&vario1, covs, constraints);
  modelPGS1.display();

  (void)vario1.dumpToNF("variopgs1.NF");
  (void)modelPGS1.dumpToNF("modelfitpgs1.NF");

  modelPGS2.fit(&vario2, covs, constraints);
  modelPGS2.display();

  (void)vario2.dumpToNF("variopgs2.NF");
  (void)modelPGS2.dumpToNF("modelfitpgs2.NF");

  RuleProp* ruleprop2;
  if (flagStationary)
    ruleprop2 = RuleProp::createFromRule(nullptr, props);
  else
    ruleprop2 = RuleProp::createFromDb(dbprop, VectorDouble());
  error = ruleprop2->fit(db, &varioparam2, 2, true);
  ruleprop2->getRule()->display();
  (void)ruleprop2->getRule()->dumpToNF("ruleFit.NF");

  modelPGS1.display();
  Vario* varioDerived = model_pgs(db, &varioparam1, ruleprop2, &modelPGS1, &modelPGS2);
  modelPGS1.display();
  varioDerived->dumpToNF("modelpgs.NF");
  varioDerived->display();

  Vario varioIndic(varioparam1);
  varioIndic.computeIndic(db, ECalcVario::VARIOGRAM);
  (void)varioIndic.dumpToNF("varioindic.NF");

  modelPGS1.display();

  delete db;
  delete dbprop;
  delete neighU;
  delete rule;
  delete ruleprop;
  delete ruleprop2;
  delete model1;
  delete model2;
  delete vario;
  delete varioDerived;

  return static_cast<int>(error);
}
