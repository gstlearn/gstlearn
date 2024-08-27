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
#include "geoslib_old_f.h"

#include "Db/Db.hpp"
#include "Db/DbHelper.hpp"
#include "Space/ASpaceObject.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 * It will be compiled but not run nor diff'ed.
 */
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("AAA_");

  defineDefaultSpace(ESpaceType::SN);

  Db* db = Db::createFillRandom(120, 2, 1);

  double mini;
  double maxi;
  double delta;
  (void)db_attribute_range(db, 2, &mini, &maxi, &delta);

  message("mini=%lf maxi=%lf delta=%lf\n", mini, maxi, delta);

  return(0);
}
