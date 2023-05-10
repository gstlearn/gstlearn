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
#include "geoslib_d.h"
#include "geoslib_f.h"

#include "Enum/ESpaceType.hpp"
#include "Enum/ECov.hpp"
#include "Enum/EKrigOpt.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  bool verbose = true;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
//  StdoutRedirect sr(sfn.str());

  // Global parameters
  defineDefaultSpace(ESpaceType::RN, 2);

  // Generate the output grid
  int nsample = 10000;
  Db* data = Db::createFillRandom(nsample);
  if (verbose) data->display();

  Timer timer;
  VectorDouble dist(nsample);
  for (int i = 0; i < nsample; i++)
    dist[i] = data->getDistance(i, 0);
  timer.displayIntervalMilliseconds("\nKriging in Unique Neighborhood");

  if (data != nullptr) delete data;

  return (0);
}
