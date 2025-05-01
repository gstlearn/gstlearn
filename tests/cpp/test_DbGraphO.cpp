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
#include "Space/ASpaceObject.hpp"

#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/DbGraphO.hpp"
#include "Basic/OptCst.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of the Db
 ** with an Oriented Graph organization
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);

  // Creating the Lines (in a 2-D space)
  defineDefaultSpace(ESpaceType::RN, 2);
  law_set_random_seed(1314);

  // tab contains the
  VectorDouble x1 = {0., 1., 2., 3., 4., 5., 6., 2., 3., 4., 5.,
                     6., 7., 3., 4., 2., 3., 4., 0., 5., 7.};
  VectorDouble x2 = {0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1.,
                     1., 1., 2., 2., 3., 3., 3., 4., 5., 4.};
  VectorDouble z1 = {1.2, 2.5, 3.6, 1.4, 0.3, 0.2, 8.2, 0.3, 3.2, 1.2, 0.4,
                     0.1, 0.3, 3.2, 4.5, 1.2, 5.2, 1.2, 1.1, 2.2, 3.3};

  VectorDouble tab;
  copy(x1.begin(), x1.end(), std::back_inserter(tab));
  copy(x2.begin(), x2.end(), std::back_inserter(tab));
  copy(z1.begin(), z1.end(), std::back_inserter(tab));

  // Constitute the connectivity graph
  NF_Triplet NF_arcs;
  law_set_random_seed(1314);
  VectorInt ind1 = {0, 1, 2, 3, 4, 5, 2, 7, 8, 9, 10, 11, 8, 13, 14, 7, 15, 16};
  VectorInt ind2 = {1,  2,  3,  4,  5,  6,  7,  8,  9,
                    10, 11, 12, 13, 14, 11, 15, 16, 17};
  int narcs = ind1.size();
  for (int iarc = 0; iarc < narcs; iarc++)
    NF_arcs.add(ind1[iarc], ind2[iarc], law_uniform(0.3, 0.9));

  DbGraphO* dbgraphO =
    DbGraphO::createFromSamples((int)x1.size(), ELoadBy::COLUMN, tab, NF_arcs,
                                {"x1", "x2", "z1"}, {"x1", "x2", "z1"});
  if (dbgraphO == nullptr) return 1;

  mestitle(1, "Reference Oriented Graph File");
   dbgraphO->display();

  // Dump into a Neutral File
  dbgraphO->dumpToNF("GraphO.ascii");

  // Read from the Neutral File and print contents again
  DbGraphO* dbgraphO2 = DbGraphO::createFromNF("GraphO.ascii");
  mestitle(1, "DbGraphO after Serialization / Deserialization");
  OptCst::define(ECst::NTROW, -1);
  OptCst::define(ECst::NTCOL, -1);
  dbgraphO2->display();

  // Check the next nodes downstream
  VectorInt iaddown = dbgraphO->getIndicesNextDown(7);
  VH::dump("Nodes next to 7 downstream", iaddown);

  // Check the next nodes upstream
  VectorInt iadup = dbgraphO->getIndicesNextUp(7);
  VH::dump("Nodes next to 7 upstream", iadup);

  // Check the end-of-stream points downwards
  VectorInt iadenddown = dbgraphO->getEndsDown();
  VH::dump("List of Node indices end-of-stream downwards", iadenddown);

  // Check the end-of-stream points upwards
  VectorInt iadendup = dbgraphO->getEndsUp();
  VH::dump("List of Node indices end-of-stream upwards", iadendup);

  // Check the orphans
  VectorInt iadorphan = dbgraphO->getOrphans();
  VH::dump("List of Node indices orphans", iadorphan);

  // Check the downstream relationship, starting from an arc number
  VectorInt order = dbgraphO->getOrderDown(3);
  VH::dump("Vertices related to #2", order);

  // Check if two nodes are connected
  message("Check if nodes 3 and 6 are connected = %d\n",
          dbgraphO->areConnected(3, 6));
  message("Check if nodes 3 and 8 are connected = %d\n",
          dbgraphO->areConnected(3, 8));

  VectorDouble cumul = dbgraphO->getCumulDown(7);
  VH::dump("Cumul of arc values starting from node 7", cumul);

  delete dbgraphO;
  delete dbgraphO2;

  return 0;
}
