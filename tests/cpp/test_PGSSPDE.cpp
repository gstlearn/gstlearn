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

#include "API/SPDE.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleProp.hpp"
#include "Model/Model.hpp"

using namespace gstlrn;

// static void _firstTest(const CovAniso* cov1,
//                        const CovAniso* cov2,
//                        std::vector<Model*>& models,
//                        Rule& rule)
// {
//   // Creating the Model(s) of the Underlying GRF(s)
//   auto* model1 = new Model();
//   model1->addCovAniso(*cov1);
//   model1->display();

//   auto* model2 = new Model();
//   model2->addCovAniso(*cov2);
//   model2->display();

//   models.push_back(model1);
//   models.push_back(model2);

//   // Creating the Rule
//   rule.resetFromNames({"S", "T", "F1", "F2", "F3"});
// }

// static void _secondTest(const CovAniso* cov1,
//                         const CovAniso* cov2,
//                         std::vector<Model*>& models,
//                         Rule& rule)
// {
//   auto* model = new Model();
//   model->addCovAniso(*cov1);
//   model->addCovAniso(*cov2);
//   model->display();

//   models.push_back(model);

//   // Creating the Rule
//   rule.resetFromNames({"S", "S", "F1", "F2", "F3"});
// }

static void _thirdTest(const CovAniso* cov1,
                       const CovAniso* cov2,
                       std::vector<Model*>& models,
                       Rule& rule)
{
  // Creating the Model(s) of the Underlying GRF(s)
  auto* model = new Model();
  model->addCovAniso(*cov1);
  model->addCovAniso(*cov2);
  model->display();

  models.push_back(model);

  // Creating the Rule
  rule.resetFromNames({"S", "S", "F1", "F2", "F3"});
}

static void _clearModels(std::vector<Model*>& models)
{
  for (size_t i = 0, n = models.size(); i < n; i++)
  {
    delete models[i];
  }
}

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix acovalistlgebra
**
*****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("test_PGSSPDE-");
  Id error  = 0;
  Id ndim   = 2;
  Id nbsimu = 3;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Prepare the Discrete process with Discretized Option
  set_test_discrete(false);

  // Prepare dimension variables
  auto ndata = 100;
  Id nx      = 101;
  double dx  = 1. / (nx - 1);

  // Prepare the output Grid
  DbGrid* grid = DbGrid::create({nx, nx}, {dx, dx});

  // Prepare facies proportions and corresponding Grid
  DbGrid* dbprop = DbGrid::create({nx, nx}, {dx, dx});
  VectorDouble props({0.2, 0.5, 0.3});
  Id nfac            = static_cast<Id>(props.size());
  VectorString names = generateMultipleNames("Props", nfac);
  for (Id ifac = 0; ifac < nfac; ifac++)
    dbprop->addColumnsByConstant(1, props[ifac], names[ifac]);
  dbprop->setLocators(names, ELoc::P, 0);

  // Prepare the input data set
  Db* dat        = Db::createFromBox(ndata, grid->getCoorMinimum(), grid->getCoorMaximum());
  VectorDouble z = VH::simulateGaussian(ndata);
  dat->addColumns(z, "Data", ELoc::Z);
  dat->display();

  // Creating the covariances involved in the Model(s) of the Underlying GRF(s)
  double range1  = 0.20;
  CovAniso* cov1 = CovAniso::createFromParam(ECov::MATERN, TEST, 1., 1., {range1, range1});

  double range2  = 0.40;
  CovAniso* cov2 = CovAniso::createFromParam(ECov::MATERN, TEST, 1., 2., {range2, range2});

  // Environment
  std::vector<Model*> models;
  Rule rule;
  RuleProp ruleprop;
  Id mode = 3;

  // IMPORTANT NOTE: the two following tests have been temporarily discraded.
  // They use the old class PGSSPDE which is now deprecated.
  // They should be reactivated as soon as the new class SPDE:
  // - allows processing a multivariate case (for PGS)
  // - allows considering Conditional PGS simulations

  // // MonoVariable PGS
  // if (mode == 0 || mode == 1)
  // {
  //   _firstTest(cov1, cov2, models, rule);
  //   ruleprop.resetFromRule(&rule, props);
  //   ruleprop.display();

  //   PGSSPDE sNonCond(models, grid, &ruleprop);
  //   law_set_random_seed(133672);
  //   sNonCond.compute(grid, 0, NamingConvention("Facies-Mono-NC"));

  //   PGSSPDE sCond(models, grid, &ruleprop, dat);
  //   law_set_random_seed(53782);
  //   sCond.compute(grid, 0, NamingConvention("Facies-Mono-CD"));

  //   _clearModels(models);
  // }

  // // Bivariate PGS
  // if (mode == 0 || mode == 2)
  // {
  //   _secondTest(cov1, cov2, models, rule);
  //   ruleprop.resetFromRule(&rule, props);
  //   ruleprop.display();

  //   PGSSPDE sNonCond(models, grid, &ruleprop);
  //   law_set_random_seed(42434);
  //   sNonCond.compute(grid, 0, NamingConvention("Facies-Multi-NC"));

  //   PGSSPDE sCond(models, grid, &ruleprop, dat);
  //   law_set_random_seed(43791);
  //   sCond.compute(grid, 0, NamingConvention("Facies-Multi-CD"));

  //   _clearModels(models);
  // }

  // Monovariate in new interface
  if (mode == 0 || mode == 3)
  {
    _thirdTest(cov1, cov2, models, rule);
    ruleprop.resetFromRule(&rule, props);
    ruleprop.display();

    (void)simPGSSPDE(dat, grid, models[0], ruleprop, nbsimu);

    _clearModels(models);
  }

  DbStringFormat dbfmt(FLAG_STATS);
  grid->display(&dbfmt);
  (void)grid->dumpToNF("pgs.NF");

  delete grid;
  delete dbprop;
  delete dat;

  return static_cast<int>(error);
}
