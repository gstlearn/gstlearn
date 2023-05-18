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
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighWork.hpp"
#include "Estimation/CalcImage.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  bool verbose = false;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
//  StdoutRedirect sr(sfn.str());

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
  NeighWork* neighW = NeighWork::create(image, neighI, image);

  Timer timer;
  krimage(image, model, neighW);
  timer.displayIntervalMilliseconds("\nKriging in Image Neighborhood");
  message("Order of magnitude of the reference implementation is 2.2Kms\n");

  if (neighI    != nullptr) delete neighI;
  if (neighW    != nullptr) delete neighW;
  if (image     != nullptr) delete image;
  if (model     != nullptr) delete model;

  return (0);
}
