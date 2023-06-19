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
//  StdoutRedirect sr(sfn.str());

  // Global parameters
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  int option = -1;

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

  // Preparing a vector of SpacePoints for the active samples in 'data'
  std::vector<SpacePoint> p1s = data->getSamplesAsSP();
  SpacePoint p2;
  VectorDouble cumul(ndat_in, 0.);
  Timer timer;

  if (option < 0 || option == 1)
  {
    // Traditional solution
    // ====================

    mestitle(1, "Traditional solution");
    timer.reset();
    for (int i = 0; i < ndat_out; i++)
    {
      dout->getSampleCoordinatesAsSP(i, p2);
      VectorDouble rhs1 = model->evalPointToDb(p2, data);
      VH::addInPlace(cumul, rhs1);
    }
    timer.displayIntervalMilliseconds("Establishing RHS", 4000);

    // Some printout for comparison
    VH::divideConstant(cumul, ndat_out);
    VH::displayRange("RHS", cumul);
  }

  if (option < 0 || option == 2)
  {
    // Semi-optimized solution
    // =======================

    mestitle(1, "Semi_optimized solution");
    VH::fill(cumul, 0.);

    timer.reset();
    for (int i = 0; i < ndat_out; i++)
    {
      dout->getSampleCoordinatesAsSP(i, p2);
      VectorDouble rhs2 = model->evalPointToDbAsSP(p2, p1s);
      VH::addInPlace(cumul, rhs2);
    }
    timer.displayIntervalMilliseconds("Establishing RHS (semi-optimized)", 1480);

    // Some printout for comparison
    VH::divideConstant(cumul, ndat_out);
    VH::displayRange("RHS", cumul);
  }

  if (option < 0 || option == 3)
  {
    // Optimized version
    // =================

    mestitle(1, "Optimized solution");
    VH::fill(cumul, 0.);

    timer.reset();
    VectorVectorDouble vecvec = model->evalCovMatrixOptim(dout, data);
    for (int i = 0; i < ndat_out; i++)
      VH::addInPlace(cumul, vecvec[i]);
    timer.displayIntervalMilliseconds("Establishing RHS (optimized)", 300);

    // Some printout for comparison
    VH::divideConstant(cumul, ndat_out);
    VH::displayRange("RHS", cumul);
  }

  // Cleaning
  if (data      != nullptr) delete data;
  if (dout      != nullptr) delete dout;
  if (model     != nullptr) delete model;

  return (0);
}
