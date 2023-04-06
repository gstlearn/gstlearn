/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_define.h"
#include "geoslib_f.h"

#include "Basic/String.hpp"
#include "Basic/File.hpp"
#include "Basic/CSVformat.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"

/**
 * This test is meant to check the CSV loading procedure
 */
int main()
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  String filename = ASerializable::getTestData("Pollution","Pollution.dat");
  Db* mydb = Db::createFromCSV(filename,CSVformat(),false);

  mydb->setLocator("X",ELoc::X,0);
  mydb->setLocator("Y",ELoc::X,1);
  mydb->setLocator("Zn",ELoc::Z);
  DbStringFormat dbfmt(FLAG_RESUME | FLAG_EXTEND | FLAG_VARS);
  mydb->display(&dbfmt);

  // Looking for duplicates
  VectorDouble dist = {0.3, 0.3};
  db_duplicate(mydb, true, dist.data());
  return 0;
}
