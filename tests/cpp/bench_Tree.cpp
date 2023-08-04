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
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Tree/Ball.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  Timer timer;
  VectorDouble vec;
  VectorDouble vecb;
  VectorDouble diff;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Tree-");

  // Global parameters
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Bench marking the Ball tree algorithm in particular
  mestitle(1, "Ball Tree Efficiency");
  int nfois = 10;
  int nech = 10000;
  VectorInt times(nfois);

  message("- Building BallTree: P(p) where p = %d * k (in ms / k)\n", nech);
  for (int ifois = 0; ifois < nfois; ifois++)
  {
    int number = nech * (ifois + 1);
    timer.reset();
    Db* data1 = Db::createFillRandom(number, ndim, 1, 0, 0., 0., VectorDouble(), VectorDouble(), VectorDouble(), 131343);
    Ball ball(data1);
    times[ifois] = timer.getIntervalMilliseconds() / (ifois + 1);

    delete data1;
  }
  VH::display("", times);

  message("- P2P: P(p1) where p1 = %d * k and P(p2=%d) (in ms / k)\n", nech, nech);

  for (int ifois = 0; ifois < nfois; ifois++)
  {
    int number = nech * (ifois + 1);
    timer.reset();
    Db* data1 = Db::createFillRandom(number, ndim, 1, 0, 0., 0., VectorDouble(), VectorDouble(), VectorDouble(), 131343);
    Db* data2 = Db::createFillRandom(nech, ndim, 1, 0, 0., 0., VectorDouble(), VectorDouble(), VectorDouble(), 413343);
    (void) migrate(data1, data2, "z", 1, VectorDouble(), true, false, true);
    times[ifois] = timer.getIntervalMilliseconds() / (ifois + 1);

    delete data1;
    delete data2;
  }
  VH::display("", times);

  message("- P2P: P(p1=%d) and P(p2) where p2 = %d * k (in ms / k)\n", nech, nech);

  for (int ifois = 0; ifois < nfois; ifois++)
  {
    int number = nech * (ifois + 1);
    timer.reset();
    Db* data1 = Db::createFillRandom(nech, ndim, 1, 0, 0., 0., VectorDouble(), VectorDouble(), VectorDouble(), 131343);
    Db* data2 = Db::createFillRandom(number, ndim, 1, 0, 0., 0., VectorDouble(), VectorDouble(), VectorDouble(), 413343);
    (void) migrate(data1, data2, "z", 1, VectorDouble(), true, false, true);
    times[ifois] = timer.getIntervalMilliseconds() / (ifois + 1);

    delete data1;
    delete data2;
  }
  VH::display("", times);

  return (0);
}
