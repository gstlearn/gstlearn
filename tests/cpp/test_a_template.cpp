/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/

#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Model/Model.hpp"
#include "Estimation/CalcGlobal.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 * It will be compiled but not run nor diff'ed.
 */
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  int ndim = 2;
  int nvar = 2;
  int ndat = 6;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Input Data File
  Db* dat = Db::createFillRandom(ndat, ndim, nvar);

  // Output Target file
  DbGrid* target = DbGrid::createFillRandom({5, 5}, 1);
  dat->setLocators({"z-1", "z-2"}, ELoc::Z);
  dat->display();
  target->setLocators({"z"}, ELoc::Z);
  target->display();

  // Model
  MatrixSymmetric* sills = MatrixSymmetric::createFromVD({2, 1, 1, 4});
  Model* model          = Model::createFromParam(ECov::EXPONENTIAL, 0.4, TEST, TEST,
                                                  VectorDouble(), *sills);
  model->setDriftIRF(0, 0);
  model->display();

  // Neighborhood (moving)
  NeighMoving* movingNeigh = NeighMoving::create(20, 100);
  NeighUnique* uniqueNeigh = NeighUnique::create();

  KrigOpt krigopt = KrigOpt();
  krigopt.setColCok({0, -1});
  OptDbg::setReference(14);

  // // In Unique Neighborhood
  // (void)kriging(dat, target, model, uniqueNeigh, true, true, false,
  //               krigopt, NamingConvention("CCKU", true, true, false));

  // In Moving Neighborhood
  (void)kriging(dat, target, model, movingNeigh, true, true, false,
                krigopt, NamingConvention("CCKM", true, true, false));

  // Various comparisons
  (void)VH::difference("CCKM.z1.estim", "z",
                       target->getColumn("CCKM.z-1.estim"),
                       target->getColumn("z"), EPSILON2, true);
  // (void)VH::difference("CCKU.z1.estim", "z",
  //                      target->getColumn("CCKU.z-1.estim"),
  //                      target->getColumn("z"), EPSILON2, true);

  // Cleaning
  delete dat;
  delete target;
  delete movingNeigh;
  delete uniqueNeigh;
  delete sills;
  delete model;
  return (0);
}
