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

/**
 * This test is meant to check the CSV loading procedure
 */
int main()
{
  String exec_dir = ASerializable::getExecDirectory();
  //String filepath("/home/fors/Documents/gstlearn/Oise/Database for kriging method Oise valley/BDD_kriggeage_Oise_elevationbottom_arcgis.csv");
  //CSVformat fmt(true, 0, ';', ',', "9999");
  //Db* mydb = new Db(filepath,true,fmt);
  //mydb->display();

  // TODO : Cross-platform way to build file path (use boost ?)
  String filepath(exec_dir + "../../../doc/data/Pollution.dat");
  Db* mydb = new Db(filepath,true,CSVformat());
  mydb->display();

  mydb->setLocator("X",ELoc::X,0);
  mydb->setLocator("Y",ELoc::X,1);
  mydb->setLocator("Zn",ELoc::Z);
  mydb->display(FLAG_RESUME | FLAG_EXTEND | FLAG_VARS);

  return 0;
}
