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

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Neigh/NeighImage.hpp"
#include "Estimation/CalcImage.hpp"

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

  // Generate the grid
  VectorInt nx = {360,240};
  DbGrid* image = DbGrid::create(nx);
  image->addColumnsByConstant(1, 1.2, "Var", ELoc::Z);
  if (verbose) image->display();

  // Create the Model
  Model* model = new Model();
  model->addCovFromParam(ECov::NUGGET, 0., 1.);
  model->addCovFromParam(ECov::SPHERICAL, 40, 2.);
  if (verbose) model->display();

  // Image Neighborhood
  NeighImage* neighI = NeighImage::create({10,10}, 3);
  if (verbose) neighI->display();

  bool flagFFT     = true;

  Timer timer;
  krimage(image, model, neighI, flagFFT);
  timer.displayIntervalMilliseconds("Kriging in Image Neighborhood", 1200);

  // Produce some stats for comparison
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_STATS, {"Filtering*"});
  image->display(dbfmt);
  delete dbfmt;

  delete neighI;
  delete image;
  delete model;

  return (0);
}
