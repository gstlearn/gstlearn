/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Db/DbStringFormat.hpp"
#include "geoslib_f.h"
#include "geoslib_define.h"
#include "Basic/String.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Basic/CSVformat.hpp"

/**
 * This test is meant to check the CSV loading procedure
 */
int main()
{
  // Standard output redirection to file
  std::stringstream sfn;
  // TODO c++17 : use #include <filesystem> to retrieve base name of __FILE__
  sfn << "test_csv" << ".out";
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
