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
#include "geoslib_define.h"

#include "Basic/File.hpp"
#include "Basic/CSVformat.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Db/DbHelper.hpp"

/**
 * This test is meant to check the CSV loading procedure
 */
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  String filename = ASerializable::getTestData("Pollution","Pollution.dat");
  Db* mydb = Db::createFromCSV(filename,CSVformat(),false);

  mydb->setLocator("X", ELoc::X, 0);
  mydb->setLocator("Y", ELoc::X, 1);
  mydb->setLocator("Zn", ELoc::Z);
  DbStringFormat dbfmt(FLAG_RESUME | FLAG_EXTEND | FLAG_VARS);
  mydb->display(&dbfmt);

  // Looking for duplicates
  VectorDouble dist = {0.3, 0.3};
  DbHelper::db_duplicate(mydb, true, dist);

  delete mydb;

  return 0;
}
