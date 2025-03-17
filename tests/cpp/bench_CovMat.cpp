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
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptCustom.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"
#include <string>

void st_bench_writing_in_matrix(int nrows, int ncols, Timer& timer)
{
  mestitle(1, "Writing into the Rectangular Covariance Matrix");
  MatrixRectangular mat(nrows, ncols);

  // consecutive writes: loop in row then col
  timer.reset();
  for (int irow = 0; irow < nrows; irow++)
    for (int icol = 0; icol < ncols; icol++)
      mat.setValue(irow, icol, 12.);
  timer.displayIntervalMilliseconds("Writing consecutively by row then by col");

  // consecutive writes: loop in col then row
  timer.reset();
  for (int icol = 0; icol < ncols; icol++)
    for (int irow = 0; irow < nrows; irow++)
      mat.setValue(irow, icol, 12.);
  timer.displayIntervalMilliseconds("Writing consecutively by col then by row");

  VectorInt rowRand = law_random_path(nrows);
  VectorInt colRand = law_random_path(ncols);

  // Writing ar random
  timer.reset();
  for (int irow = 0; irow < nrows; irow++)
    for (int icol = 0; icol < ncols; icol++)
      mat.setValue(rowRand[irow], colRand[icol], 12.);
  timer.displayIntervalMilliseconds("Writing randomly by row then by col");

  // Writing ar random
  timer.reset();
  for (int icol = 0; icol < ncols; icol++)
    for (int irow = 0; irow < nrows; irow++)
      mat.setValue(rowRand[irow], colRand[icol], 12.);
  timer.displayIntervalMilliseconds("Writing randomly by col then by row");
}

  /****************************************************************************/
  /*!
   ** Main Program
   **
   *****************************************************************************/
  int main(int argc, char* argv[])
  {
    bool verbose = false;

    std::stringstream sfn;
    sfn << gslBaseName(__FILE__) << ".out";
    StdoutRedirect sr(sfn.str(), argc, argv);

    // Global parameters
    int mode = 3;
    int ndim = 2;
    defineDefaultSpace(ESpaceType::RN, ndim);

    // Generate the input data base
    int nall = 100;
    Db* dbin = Db::createFillRandom(nall, ndim);
    dbin->addSelectionRandom(0.9);
    int ndat = dbin->getNSample(true);
    if (verbose) dbin->display();

    // Generate the output data base
    int nout  = 100000;
    Db* dbout = Db::createFillRandom(nout, ndim);
    if (verbose) dbout->display();

    // Create the Model
    double range = 0.6;
    double sill  = 1.2;
    Model* model = Model::createFromParam(ECov::SPHERICAL, range, sill);
    if (verbose) model->display();

    // Printout
    message("RHS between:\n");
    message("- each active sample (%d out of %d) of the input data base\n", ndat, nall);
    message("- each one of the %d target sites\n", nout);
    message("(For checking purpose, a Selection has been added)\n");
    message("Statistics are provided on the averaged RHS\n");

    // Core allocation of common variables
    SpacePoint p2(model->getSpace());
    VectorDouble cumul(ndat, 0.);
    VectorDouble rhs1;
    Timer timer;

    if (mode == 0 || mode == 1)
    {
      // Traditional solution
      // ====================

      mestitle(1, "Traditional solution");
      message("Double loop on the input and output points\n");
      timer.reset();
      for (int i = 0; i < nout; i++)
      {
        dbout->getSampleAsSPInPlace(p2, i);
        model->evalPointToDb(rhs1, p2, dbin);
        VH::addInPlace(cumul, rhs1);
      }
      timer.displayIntervalMilliseconds("Establishing RHS", 3900);

      // Some printout for comparison
      VH::divideConstant(cumul, nout);
      VH::dumpRange("", cumul);
    }

    if (mode == 0 || mode == 2)
    {
      // Semi-optimized solution
      // =======================

      mestitle(1, "Semi_optimized solution");
      message("Input samples are pre-transformed into vector of space points\n");
      message("Simple loop between each target and the previous vector\n");
      VH::fill(cumul, 0.);

      timer.reset();

      // Preparing a vector of SpacePoints for the active samples in 'data'
      // for this usage, the list of SP can be reduced to the active samples only
      std::vector<SpacePoint> p1s;
      dbin->getSamplesAsSP(p1s, model->getSpace(), true);
      VectorDouble rhs2;

      for (int i = 0; i < nout; i++)
      {
        dbout->getSampleAsSPInPlace(p2, i);
        model->evalPointToDbAsSP(rhs2, p1s, p2);
        VH::addInPlace(cumul, rhs2);
      }
      timer.displayIntervalMilliseconds("Establishing RHS (semi-optimized)", 600);

      // Some printout for comparison
      VH::divideConstant(cumul, nout);
      VH::dumpRange("", cumul);
    }

    if (mode == 0 || mode == 3)
    {
      // Optimized version
      // =================

      mestitle(1, "Optimized solution");
      message("Input samples are pre-transformed into vector of (anisotropic) space points\n");
      message("Simple loop between each target and the previous vector\n");
      model->setOptimEnabled(true);

      timer.reset();
      OptCustom::define("OptimCovMat", mode);
      MatrixRectangular mat = model->evalCovMat(dbin, dbout);
      timer.displayIntervalMilliseconds("Establishing RHS V" + std::to_string(int(mode)));

      // Some printout for comparison
      VH::fill(cumul, 0.);
      for (int i = 0; i < nout; i++)
        VH::addInPlace(cumul, mat.getColumn(i));
      VH::divideConstant(cumul, nout);
      VH::dumpRange("", cumul);
    }

    if (mode == 0 || mode == 4)
    {
      // Measure the difference between consecutive access vs. random access
      st_bench_writing_in_matrix(nall, nout, timer);
    }

    // Cleaning
    delete dbin;
    delete dbout;
    delete model;

    return (0);
  }
