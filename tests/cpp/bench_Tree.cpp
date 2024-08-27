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
#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Tree/Ball.hpp"

/****************************************************************************/
/*!
 ** Main Program
 ** This is meant to test the time improvement using BallTree search or not.
 ** It is illustrated within the 'ligrate' algorithm.
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
  bool flag_stats = false; // This is set to FALSE in order to avoid Diff in Time

  // Bench marking the Ball tree algorithm in particular
  mestitle(1, "Ball Tree Efficiency");
  int nfois = 10;
  int nech = 10000;
  VectorInt times(nfois);
  if (! flag_stats)
    message("To get statistics on Ball Tree Efficiency, turn 'flag_stats' to TRUE\n");

  message("- Building BallTree: Db(n = %d * k) (in ms per k)\n", nech);
  for (int ifois = 0; ifois < nfois; ifois++)
  {
    int number = nech * (ifois + 1);
    timer.reset();
    Db *data1 = Db::createFillRandom(number, ndim, 1, 0, 0, 0., 0.,
                                     VectorDouble(), VectorDouble(),
                                     VectorDouble(), 131343);
    Ball ball(data1);
    times[ifois] = timer.getIntervalMilliseconds() / (ifois + 1);

    delete data1;
  }
  if (flag_stats) VH::display("", times);

  message("- Migrate P2P: from Db1(n = %d * k) to Db2(n = %d) (in ms per k)\n", nech, nech);

  for (int ifois = 0; ifois < nfois; ifois++)
  {
    int number = nech * (ifois + 1);
    timer.reset();
    Db *data1 = Db::createFillRandom(number, ndim, 1, 0, 0, 0., 0.,
                                     VectorDouble(), VectorDouble(),
                                     VectorDouble(), 131343);
    Db *data2 = Db::createFillRandom(nech, ndim, 1, 0, 0, 0., 0.,
                                     VectorDouble(), VectorDouble(),
                                     VectorDouble(), 413343);
    (void) migrate(data1, data2, "z", 1, VectorDouble(), true, false, true);
    times[ifois] = timer.getIntervalMilliseconds() / (ifois + 1);

    delete data1;
    delete data2;
  }
  if (flag_stats) VH::display("", times);

  message("- Migrate P2P: from Db1(n = %d) to Db2(n = %d * k) (in ms per k)\n", nech, nech);

  for (int ifois = 0; ifois < nfois; ifois++)
  {
    int number = nech * (ifois + 1);
    timer.reset();
    Db *data1 = Db::createFillRandom(nech, ndim, 1, 0, 0, 0., 0.,
                                     VectorDouble(), VectorDouble(), VectorDouble(), 131343);
    Db *data2 = Db::createFillRandom(number, ndim, 1, 0, 0, 0., 0.,
                                     VectorDouble(), VectorDouble(), VectorDouble(), 413343);
    (void) migrate(data1, data2, "z", 1, VectorDouble(), true, false, true);
    times[ifois] = timer.getIntervalMilliseconds() / (ifois + 1);

    delete data1;
    delete data2;
  }
  if (flag_stats) VH::display("", times);

  return (0);
}
