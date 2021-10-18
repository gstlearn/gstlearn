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

/**
 * This test is meant to check the CSV loading procedure
 */
int main()
{
  String filepath("doc/data/Pollution.dat");
  Db* mydb = new Db(filepath,true,CSVformat());
  mydb->setLocator(VectorString({"X","Y"}),ELoc::X);
  mydb->setLocator("Zn",ELoc::Z);
  mydb->display(FLAG_RESUME | FLAG_EXTEND | FLAG_VARS);

  return 0;
}
