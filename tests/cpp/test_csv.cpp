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
#include "Db/Db.hpp"
#include "Basic/CSVformat.hpp"


  // TODO : Cross-platform way to build file path (use boost ?)
String getTestData(const String& filename)
{
  String exec_dir = ASerializable::getExecDirectory();
  // This path is compatible with CMake generation
  String filepath(exec_dir + "../../doc/data/" + filename);

  return filepath;
}



/**
 * This test is meant to check the CSV loading procedure
 */
int main()
{
  String filepath = getTestData("Pollution.dat");
  Db* mydb = new Db(filepath,true,CSVformat());

  mydb->setLocator("X",ELoc::X,0);
  mydb->setLocator("Y",ELoc::X,1);
  mydb->setLocator("Zn",ELoc::Z);
  DbStringFormat dbfmt;
  dbfmt.setParams(FLAG_RESUME | FLAG_EXTEND | FLAG_VARS);
  mydb->display(&dbfmt);

  // Looking for duplicates
  VectorDouble dist = {0.3, 0.3};
  db_duplicate(mydb, true, dist.data());
  return 0;
}
