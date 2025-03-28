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
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Db/RankHandler.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the facility offered by RankHandler class
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("TestRkHdl-");
  int seed = 10355;
  law_set_random_seed(seed);
  OptCst::define(ECst::NTROW, -1);
  OptCst::define(ECst::NTCOL, -1);

  // Creating the Db
  int ndat = 20;
  int ndim = 2;
  int nvar = 3;
  double selRatio = 0.2;
  VectorDouble heteroRatio = {0.3, 0.2, 0.1};
  Db* db = Db::createFillRandom(ndat, ndim, nvar, 0, 0, 0., selRatio, heteroRatio);
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_RESUME | FLAG_ARRAY);
  db->display(dbfmt);

  // Instanciate the RankHandler
  RankHandler* rkhd = new RankHandler(db, true, true, true, true);

  // Creating the vector of elligible samples
  // It is constructed so as to involve:
  // - masked samples
  // - heterotopic sample
  // Warning: Db contents should not be mofiied anymore.
  VectorInt nbgh = {0, 8, 10, 12, 13};
  VH::dump("List of ranks for Elligible samples", nbgh);

  // Define the Rank Handler for the previous list
  rkhd->defineSampleRanks(nbgh);

  // Using some features of the RankHandler class
  message("Number of elligible samples: %d\n", rkhd->getNumber());
  message("Total count of sample references = %d\n", rkhd->getTotalCount());
  for (int ivar = 0; ivar < nvar; ivar++)
    message("Number of references for Variable %d = %d\n", ivar, rkhd->getCount(ivar));

  // Complete dump
  rkhd->dump(true);

  delete db;
  delete dbfmt;
  delete rkhd;

  return 0;
}

