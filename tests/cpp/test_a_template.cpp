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
#include "Variogram/Variograms.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  auto* db = Db::createFillRandom(10);
  auto* vcloud = vcloudFromDb(db);
  vcloud->getStatsAsTable().display();
  delete vcloud;

  // Define an Interval
  auto distmax = db->getExtensionDiagonal() / 3.;
  auto varmin = db->getVariance(db->getNameByLocator(ELoc::Z));
  vcloud = vcloudFromDb(
    db, nullptr, TEST, ECalcVario::VARIOGRAM, true, TEST, 100, 100, nullptr,
    distmax, varmin);
  vcloud->getStatsAsTable().display();
  delete vcloud;
  return 0;
}
