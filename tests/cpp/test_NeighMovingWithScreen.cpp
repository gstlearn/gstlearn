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
#include <cassert>

#include "Db/Db.hpp"
#include "Neigh/NeighMovingWithScreens.hpp"
#include "Neigh/NeighWork.hpp"
#include "Basic/Vector.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program demonstrates the capabilities of NeighMovingWithScreens class
 **
 *****************************************************************************/


int main(int /*argc*/, char */*argv*/[])
{
  auto *dbin = Db::createFromSamples(
    4, ELoadBy::SAMPLE, 
    {
      0, 0, 0,  // P0
      1, 0, 0,  // P1
      0, 1, 1,  // P2
      1, 1, 1   // P3
    },
    {"x1", "x2", "z"}, {"x1", "x2", "z1"}, 0
  );
  auto *dbout = Db::createFromSamples(
    4, ELoadBy::SAMPLE, 
    {
      0, 0.5,  // Q0
      0.5, 0,  // Q1
      1, 0.5,  // Q2
      0.5, 1   // Q3
    },
    {"x1", "x2"}, {"x1", "x2"}, 0
  );
  const auto screens = VectorDouble{
      0.2, 0.5, 0.8, 0.5,  // F1: (0.2, 0.5) - (0.8, 0.5)
      0.2, 0.4, 0.2, 0.6,  // F2: (0.2, 0.4) - (0.2, 0.6)
      0.8, 0.4, 0.8, 0.6,  // F3: (0.8, 0.4) - (0.8, 0.6)
      0.1, 0.1, 0.1, 0.1   // Invalid, should not be added 
  };
  // F1 intersects (P0, Q3), (P1, Q3), (P2, Q1), (P3, Q1)
  // F2 intersects (P1, Q0), (P3, Q0) (Note: intersect at F2 vertices)
  // F3 intersects (P0, Q2), (P2, Q2) (Note: intersect at F3 vertices)
  auto *neigh = NeighMovingWithScreens::create(2, false, screens, 24);
  assert(neigh != nullptr);
  assert(neigh->getType() == ENeigh::MOVING_WITH_SCREENS);
  assert(neigh->getScreens().size() == 3);

  NeighWork nwork{dbin, neigh};
  VectorInt ranks, ranks_ref;
  // Q0
  ranks = nwork.select(dbout, 0);
  ranks_ref = VectorInt{0, 2};
  if (!(ranks == ranks_ref)) {throw std::exception();}
  // Q1
  ranks = nwork.select(dbout, 1);
  ranks_ref = VectorInt{0, 1};
  if (!(ranks == ranks_ref)) {throw std::exception();}
  // Q2
  ranks = nwork.select(dbout, 2);
  ranks_ref = VectorInt{1, 3};
  if (!(ranks == ranks_ref)) {throw std::exception();}
  // Q3
  ranks = nwork.select(dbout, 3);
  ranks_ref = VectorInt{2, 3};
  if (!(ranks == ranks_ref)) {throw std::exception();}

  delete dbin;
  delete dbout;
  delete neigh;
  return 0;
}
