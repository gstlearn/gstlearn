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
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Space/SpaceRN.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  DECLARE_UNUSED(argc);
  DECLARE_UNUSED(argv);

  // Define the two Dbs in different space dimensions
  Id ndat   = 100;
  Id ndim2  = 2;
  Id ndim3  = 3;
  auto* db2 = Db::createFillRandom(ndat, ndim2, 1);
  auto* db3 = Db::createFillRandom(ndat, ndim3, 1);

  // Define the Moving neighborhood
  Id nmaxi      = 10;
  double radius = 0.5;
  auto* neigh   = NeighMoving::create(false, nmaxi, radius);

  // Attach the neighborhood to the first Db (Space Dimension 2)
  neigh->attach(db2, db2);
  neigh->display();

  VectorInt nbgh;
  std::vector<SpacePoint> pvec;

  // Perform a neighborhood selection on Db2
  Id itarget2 = 4;
  neigh->select(itarget2, nbgh);
  VH::dump("For Target Point", nbgh);
  auto space2 = SpaceRN::create(ndim2);
  db2->getSamplesFromNbghAsSP(pvec, space2, nbgh);
  for (Id i = 0; i < nmaxi; i++)
    pvec[i].display();

  // Attach the neighborhood to the second Db (Space Dimension 3)
  neigh->attach(db3, db3);
  neigh->display();

  // Perform a neighborhood selection on Db3
  Id itarget3 = 34;
  neigh->select(itarget3, nbgh);
  VH::dump("For Target Point", nbgh);
  auto space3 = SpaceRN::create(ndim3);
  db3->getSamplesFromNbghAsSP(pvec, space3, nbgh);
  for (Id i = 0; i < nmaxi; i++)
    pvec[i].display();

  delete db2;
  delete db3;
  return 0;
}
