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
#include "Enum/ESpaceType.hpp"
#include "Enum/ECov.hpp"
#include "Enum/EKrigOpt.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  bool verbose = false;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Global parameters
  defineDefaultSpace(ESpaceType::RN, 2);

  // Generate the output grid
  int nsample = 1000000;
  Db* data = Db::createFillRandom(nsample);
  if (verbose) data->display();

  Timer timer;
  VectorDouble dist(nsample);
  for (int i = 0; i < nsample; i++)
    dist[i] = data->getDistance(i, 0);
  timer.displayIntervalMilliseconds("Kriging in Unique Neighborhood", 310);

  // Produce some statistics for comparison
  DbStringFormat *dbfmt = DbStringFormat::create(FLAG_STATS, {"z"});
  data->display(dbfmt);

  delete dbfmt;

  if (data != nullptr) delete data;

  return (0);
}
