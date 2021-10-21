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
#include "Db/Db.hpp"
#include "Basic/CSVformat.hpp"
#include "Utility.hpp"
#include "geoslib_f.h"

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
  mydb->displayMore(FLAG_RESUME | FLAG_EXTEND | FLAG_VARS);

  // Looking for duplicates
  VectorDouble dist = {0.3, 0.3};
  db_duplicate(mydb, true, dist.data());
  return 0;
}
