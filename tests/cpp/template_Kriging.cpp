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

// This case study is meant to demonstrate a Kriging procedure:
// - either by using the corresponding API
// - or by solving the Kriging problem by hand

#include "geoslib_f.h"

#include "Enum/ESpaceType.hpp"

#include "Basic/OptCst.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Estimation/KrigOpt.hpp"
#include "Estimation/KrigingAlgebra.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Global parameters
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  OptCst::define(ECst::NTROW, -1);
  OptCst::define(ECst::NTCOL, -1);
  bool verbose = true;

  // Generate the input data (in the 1x1 square) with 'nvar' variables and 'nfex' external drifts
  int ndat = 10;
  int nvar   = 2;
  int nfex   = 2;
  int seedin = 13227;
  Db* dbin   = Db::createFillRandom(ndat, ndim, nvar, nfex, 0, 0., 0.,
                                    VectorDouble(), VectorDouble(), VectorDouble(), seedin);
  if (verbose) dbin->display();

  // Generate the output data set
  int nout    = 20;
  int seedout = 134484;
  Db* dbout   = Db::createFillRandom(nout, ndim, 0, nfex, 0, 0., 0.,
                                     VectorDouble(), VectorDouble(), VectorDouble(), seedout);
  if (verbose) dbout->display();

  // Create the Model
  int order               = 0;
  std::vector<ECov> types = {ECov::NUGGET, ECov::SPHERICAL};
  Model* model            = Model::createFillRandom(ndim, nvar, types, 1., order, nfex);
  if (verbose) model->display();

  // Creating a Moving Neighborhood
  NeighUnique* neigh = NeighUnique::create(false);
  if (verbose) neigh->display();

  // Perform Kriging using the API
  mestitle(0, "Kriging using API");
  (void)kriging(dbin, dbout, model, neigh, EKrigOpt::POINT, true, true, true);
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_ARRAY, {"Kriging*"});
  dbout->display(dbfmt); // Always print the results at all target sites

  // Perform Kriging 'by hand'
  mestitle(0, "Kriging constructed 'by hand'");
  MatrixSquareSymmetric Sigma00 = model->eval0Mat();
  if (verbose) Sigma00.dumpStatistics("C00 Matrix");
  MatrixSquareSymmetric Sigma   = model->evalCovMatSym(dbin);
  if (verbose) Sigma.dumpStatistics("LHS: Covariance part");
  MatrixRectangular X           = model->evalDriftMat(dbin);
  if (verbose) X.dumpStatistics("LHS: Drift part");

  VectorVectorInt sampleRanks   = dbin->getSampleRanks();
  VectorDouble Z                = dbin->getValuesByRanks(sampleRanks);

  KrigingAlgebra Kcalc;
  Kcalc.setData(&Z, &sampleRanks);
  Kcalc.setLHS(&Sigma, &X);
  Kcalc.setVariance(&Sigma00);
  MatrixRectangular Sigma0;
  MatrixRectangular X0;
  KrigOpt krigopt;

  // Loop on the target sites
  VectorDouble result;
  for (int iout = 0; iout < nout; iout++)
  {
    if (model->evalCovMatRHSInPlaceFromIdx(Sigma0, dbin, dbout, sampleRanks, iout, krigopt, false)) break;
    if (verbose && iout == 0) Sigma0.dumpStatistics("RHS(target:1): Covariance part");  
    if (model->evalDriftMatByTarget(X0, dbout, iout, krigopt)) break;
    if (verbose && iout == 0) X0.dumpStatistics("RHS(target:1): Drift part");

    Kcalc.setRHS(&Sigma0, &X0);
    result.clear();
    VH::concatenateInPlace(result, Kcalc.getEstimation());
    VH::concatenateInPlace(result, Kcalc.getStdv());
    VH::concatenateInPlace(result, Kcalc.getVarianceZstar()); 
    VH::dump("Sample " + std::to_string(iout+1), result, false); // Print results at all target sites
  }

  // Free classes
  delete dbin;
  delete dbout;
  delete model;
  delete neigh;
  delete dbfmt;

  return (0);
}
