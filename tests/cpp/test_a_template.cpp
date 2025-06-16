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
  defineDefaultSpace(ESpaceType::RN, ndim);

  Db* dat = Db::createFillRandom(20, 2, 2);
  dat->display();

  DbGrid* target = DbGrid::createFillRandom({10,10}, 1);
  target->display();
  dat->setLocators({"z-1", "z-2"}, ELoc::Z);
  target->setLocators({"z"}, ELoc::Z);

  VarioParam* varioparam = VarioParam::createMultiple(2, 30, 10);
  Vario vario_raw2dir    = Vario(*varioparam);
  (void) vario_raw2dir.compute(dat);

  Model fitmod_raw = Model();
  (void) fitmod_raw.fit(&vario_raw2dir,
                        {ECov::NUGGET, ECov::EXPONENTIAL, ECov::CUBIC, ECov::LINEAR});
  fitmod_raw.setDriftIRF(0, 0);
  fitmod_raw.display();

  KrigOpt krigopt          = KrigOpt();
  NeighUnique* uniqueNeigh = NeighUnique::create();
  NeighMoving* movingNeigh = NeighMoving::create(20, 100);

  krigopt.setColCok({1, -1});

  // this one throws an error
  (void)kriging(dat, target, &fitmod_raw, uniqueNeigh,true,false,false,
                krigopt, NamingConvention("COLCOK"));

  // this one runs but results are not correct
  (void)kriging(dat, target, &fitmod_raw, movingNeigh,true,false,false,
                krigopt, NamingConvention("COLCOK"));
  return (0);
}
