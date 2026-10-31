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
#include "Basic/ASerializable.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Estimation/Estimations.hpp"
#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"
#include "geoslib_define.h"
#include "utils.hpp"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Unless you test a specific feature, bring it to its minimal expansion
  // e.g. the few lines below
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-");

  auto fileNF = getTestData("Scotland", "Scotland_Temperatures.NF");
  auto* dat = Db::createFromNF(fileNF);
  dat->setLocators(VectorString{"January_temp", "Elevation"}, ELoc::Z);

  fileNF = getTestData("Scotland", "Scotland_Elevations.NF");
  auto* target = DbGrid::createFromNF(fileNF);

  auto* varioparam = VarioParam::createMultiple(2, 30, 10);
  auto* varioexp2var = Vario::create(*varioparam);
  (void)varioexp2var->compute(dat);

  Model fitmod2var = Model();
  (void)fitmod2var.fit(
    varioexp2var, {ECov::NUGGET, ECov::EXPONENTIAL, ECov::CUBIC, ECov::LINEAR});
  (void)fitmod2var.setDriftIRF(0, 0);

  dat->display();
  dat->getStatsAsTable().display();
  (void)kriging(dat, target, &fitmod2var);
  target->getStatsAsTable().display();
  return 0;
}
