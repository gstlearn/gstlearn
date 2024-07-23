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
int main(int argc, char *argv[])
{
  bool verbose = false;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Global parameters
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  int option = -1;

  // Generate the input data base
  int nall = 100;
  Db* dbin = Db::createFillRandom(nall, ndim);
  dbin->addSelectionRandom(0.9);
  int ndat = dbin->getSampleNumber(true);
  if (verbose) dbin->display();

  // Generate the output data base
  int nout = 100000;
  Db* dbout = Db::createFillRandom(nout, ndim);
  if (verbose) dbout->display();

  // Create the Model
  double range = 0.6;
  double sill = 1.2;
  Model* model = Model::createFromParam(ECov::SPHERICAL, range, sill);
  if (verbose) model->display();

  // Printout
  message("RHS between:\n");
  message("- each one of the %d target sites\n",nout);
  message("- all (active) samples (%d) of the input data base\n",ndat);
  message("(For checking purpose, a Selection has been added)\n");
  message("Statistics are provided on the averaged RHS\n");

  // Preparing a vector of SpacePoints for the active samples in 'data'
  // for this usage, the list of SP can be reduced to the active samples only
  std::vector<SpacePoint> p1s = dbin->getSamplesAsSP(true);
  SpacePoint p2;
  VectorDouble cumul(ndat, 0.);
  Timer timer;

  if (option < 0 || option == 1)
  {
    // Traditional solution
    // ====================

    mestitle(1, "Traditional solution");
    message("Double loop on the input and output points\n");
    timer.reset();
    for (int i = 0; i < nout; i++)
    {
      dbout->getSampleAsSPInPlace(i, p2);
      VectorDouble rhs1 = model->evalPointToDb(p2, dbin);
      VH::addInPlace(cumul, rhs1);
    }
    timer.displayIntervalMilliseconds("Establishing RHS", 3900);

    // Some printout for comparison
    VH::divideConstant(cumul, nout);
    VH::displayRange("", cumul);
  }

  if (option < 0 || option == 2)
  {
    // Semi-optimized solution
    // =======================

    mestitle(1, "Semi_optimized solution");
    message("Input samples are pre-transformed into vector of space points\n");
    message("Simple loop between each target and the previous vector\n");
    VH::fill(cumul, 0.);

    timer.reset();
    for (int i = 0; i < nout; i++)
    {
      dbout->getSampleAsSPInPlace(i, p2);
      VectorDouble rhs2 = model->evalPointToDbAsSP(p1s, p2);
      VH::addInPlace(cumul, rhs2);
    }
    timer.displayIntervalMilliseconds("Establishing RHS (semi-optimized)", 600);

    // Some printout for comparison
    VH::divideConstant(cumul, nout);
    VH::displayRange("", cumul);
  }

  if (option < 0 || option == 3)
  {
    // Optimized version
    // =================

    mestitle(1, "Optimized solution");
    message("Input samples are pre-transformed into vector of (anisotropic) space points\n");
    message("Simple loop between each target and the previous vector\n");
    VH::fill(cumul, 0.);

    timer.reset();
    MatrixRectangular matvec = model->evalCovMatrixOptim(dbin, dbout);
    for (int i = 0; i < nout; i++)
      VH::addInPlace(cumul, matvec.getColumn(i));
    timer.displayIntervalMilliseconds("Establishing RHS (optimized)", 300);

    // Some printout for comparison
    VH::divideConstant(cumul, nout);
    VH::displayRange("", cumul);
  }

  // Cleaning
  if (dbin      != nullptr) delete dbin;
  if (dbout     != nullptr) delete dbout;
  if (model     != nullptr) delete model;

  return (0);
}
