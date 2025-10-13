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
#include "Enum/EFormatNF.hpp"
#include "Enum/ESpaceType.hpp"
#include "Space/ASpaceObject.hpp"

#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Db/Db.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the use of Neutral files
 ** This is performed on a standard Db for illustration
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("bench_NF-");

  // Creating the Lines (in a 2-D space)
  defineDefaultSpace(ESpaceType::RN, 2);
  Timer timer;

  // Creating the Db
  Id ndat  = 1000000;
  Id ndim  = 2;
  Id nvar  = 10;
  Db* db1  = Db::createFillRandom(ndat, ndim, nvar);
  db1->display();

  // Saving as a Neutral file
  timer.reset();
  db1->dumpToNF("Db.NF.ascii", EFormatNF::ASCII);
  timer.displayIntervalMilliseconds("Dump into NF");

  // Reading the neutral file
  timer.reset();
  Db* db2 = Db::createFromNF("Db.NF.ascii");
  db2->display();
  timer.displayIntervalMilliseconds("Load from NF");

  delete db1;
  delete db2;

  return 0;
}
