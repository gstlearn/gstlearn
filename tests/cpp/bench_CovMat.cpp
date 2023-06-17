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
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"

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
  StdoutRedirect sr(sfn.str());

  // Global parameters
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Generate the input data base
  int ndat_in = 100;
  Db* data = Db::createFillRandom(ndat_in, ndim);
  if (verbose) data->display();

  // Generate the output data base
  int ndat_out = 100000;
  Db* dout = Db::createFillRandom(ndat_out, ndim);
  if (verbose) dout->display();

  // Create the Model
  double range = 0.6;
  double sill = 1.2;
  Model* model = Model::createFromParam(ECov::SPHERICAL, range, sill);
  if (verbose) model->display();

  // Printout
  message("RHS between:\n");
  message("- each one of the %d target sites\n",ndat_out);
  message("- all samples (%d) of the input data base\n",ndat_in);
  message("Statistics are provided on the averaged RHS\n");

  // Traditional solution
  SpacePoint p1;
  VectorDouble cumul(ndat_in, 0.);
  Timer timer;
  for (int i = 0; i < ndat_out; i++)
  {
    dout->getSampleCoordinatesAsSP(i, p1);
    VectorDouble vec = model->evalPointToDb(p1, data);
    VH::addInPlace(cumul, vec);
  }
  timer.displayIntervalMilliseconds("Establishing RHS", 3800);

  // Some printout for comparison
  VH::divideConstant(cumul, ndat_out);

  VH::displayRange("RHS", cumul);

  // Convert the contents of the data base as a vector of Space Points
  std::vector<SpacePoint> pvec = data->getSamplesAsSP();

  // Checking using the optimized version
  VH::fill(cumul,  0.);

  timer.reset();
  model->getCovAnisoList()->preProcess(pvec);
  VectorVectorDouble vecvec = model->evalCovMatrixOptim(dout, data);
  for (int i = 0; i < ndat_out; i++)
    VH::addInPlace(cumul, vecvec[i]);
  model->getCovAnisoList()->cleanPreProcessInfo();
  timer.displayIntervalMilliseconds("Establishing RHS (optimized)", 3800);

  // Some printout for comparison
  VH::divideConstant(cumul, ndat_out);
  VH::displayRange("RHS", cumul);

  if (data      != nullptr) delete data;
  if (dout      != nullptr) delete dout;
  if (model     != nullptr) delete model;

  return (0);
}
