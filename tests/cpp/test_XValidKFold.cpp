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
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/VectorNumT.hpp"
#include "Db/Db.hpp"
#include "Enum/ECst.hpp"
#include "Enum/ELoadBy.hpp"
#include "Enum/ESpaceType.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Space/ASpaceObject.hpp"
#include "Tree/Ball.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 ** This is meant to test XValidation neighborhood with K-Fold option
 ** using Ball tree
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("test_XValidKFold-");

  // Global parameters
  defineDefaultSpace(ESpaceType::RN, 2);
  OptCst::define(ECst::NTCAR, 3);

  // Main parameters
  VectorString names({"X", "Y", "S", "C", "V"});
  VectorDouble vals({0.5, 1.0, 1.5, 1.0, 1.5, 8.0, 9.0,   // X
                     0.5, 1.0, 1.0, 1.5, 1.5, 8.0, 9.0,   // Y
                     1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,   // Sel
                     1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0,   // Code
                     0.1, 0.2, 0.3, 0.4, 0.5, 1.6, 1.7}); // Variable

  // Number of data
  Id ndat = 7;

  // Neighborhood parameters
  Id nmaxi      = 2;
  Id leaf_size  = 5;
  double radius = 1;

  // Generate the input data base
  Db* data = Db::createFromSamples(ndat, ELoadBy::COLUMN, vals, names);
  data->setLocators({"X", "Y"}, ELoc::X);
  data->setLocator("S", ELoc::SEL);
  data->setLocator("C", ELoc::C);
  data->setLocator("V", ELoc::Z);
  data->display();

  // Moving Neighborhood
  NeighMoving* neigh = NeighMoving::create(true, nmaxi, radius);
  neigh->setFlagXvalid(true);
  neigh->setFlagKFold(true);
  neigh->setBallSearch(true, leaf_size);
  neigh->attach(data, data); // XValidation (dbout = dbin)

  // Display the ranks
  VectorInt ranks;
  mestitle(0, "Neigh ranks for each target:");
  for (Id i = 0; i < ndat; i++)
  {
    neigh->select(i, ranks);
    String str = toStrFormat("Neigh ranks [%d]:", i);
    if (ranks.empty())
      message("%s None\n", str.c_str());
    else
      printVector(ranks, str);
  }

  delete neigh;
  delete data;
  return (0);
}
