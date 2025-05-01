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

// This test is mean to check the Factorial Kriging Analysis On Grid

#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCst.hpp"
#include "Neigh/NeighImage.hpp"
#include "Estimation/CalcImage.hpp"

int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Global parameters
  int ndim = 2;
  law_set_random_seed(32131);
  defineDefaultSpace(ESpaceType::RN, ndim);
  DbStringFormat* dbfmt;
  OptCst::define(ECst::NTROW, -1);

  // Parameters
  bool debug         = false;
  bool verbose       = true;
  int nx             = 100;
  int ny             = 100;
  int nvar           = 3;
  int skip           = 0;
  bool flagSK        = true;
  bool flagFFT       = true;
  VectorInt radius   = {10, 10};
  VectorDouble means = {0., 100., -50.}; // Must be adapted to the number of variables

  // Generate the target file
  DbGrid* db = DbGrid::createFillRandom({nx, ny}, nvar, 0, 0, 0., 0., VectorDouble(), means);
  db->display();

  // Create the Model
  int order    = (flagSK) ? -1 : 0;
  Model* model = Model::createFillRandom(ndim, nvar,
                                         {ECov::NUGGET, ECov::SPHERICAL},
                                         1., order);
  model->setCovFiltered(0, true);
  if (flagSK) model->setMeans(means);
  model->display();

  // Neighborhood
  NeighImage* neigh = NeighImage::create(radius, skip);
  neigh->display();

  // Define the verbose option
  if (debug) OptDbg::setReference(1);

  // Test on Collocated CoKriging in Unique Neighborhood
  (void) krimage(db, model, neigh, flagFFT, verbose);
  dbfmt = DbStringFormat::create(FLAG_STATS, {"Filtering*"});
  db->display(dbfmt);

  // Free pointers

  delete neigh;
  delete db;
  delete dbfmt;
  delete model;

  return (0);
}
