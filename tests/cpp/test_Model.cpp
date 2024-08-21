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
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Drifts/DriftM.hpp"
#include "Basic/File.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/OptCst.hpp"
#include "Enum/ESpaceType.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Model-");
  int seed = 10355;
  int ndim = 2;
  int nvar = 1;
  law_set_random_seed(seed);

  ///////////////////////
  // Creating the Db
  auto nx = { 2, 3 };
  auto x0 = { 5.2, 8.3 };
  auto dx = { 1.3, 0.6 };
  DbGrid* workingDbc = DbGrid::create(nx, dx, x0);

  // Building the Covariance Context
  CovContext ctxt(nvar, ndim);

  ///////////////////////
  // Creating the Model
  Model modellmc = Model(ctxt);
  // Build the List of Covariances
  CovLMC covlmc = CovLMC(ctxt.getSpace());
  // Build the Elementary Covariances
  CovAniso cov1 = CovAniso(ECov::CUBIC,ctxt);
  cov1.setRanges({1.2,2.1});
  cov1.setSill(1.5);
  covlmc.addCov(&cov1);
  CovAniso cov2 = CovAniso(ECov::NUGGET,ctxt);
  cov2.setSill(0.5);
  covlmc.addCov(&cov2);
  // Assembling the Model
  modellmc.setCovList(&covlmc);
  modellmc.display();

  // Building the Covariance Matrix
  MatrixSquareSymmetric result = modellmc.evalCovMatrixSymmetric(workingDbc);
  result.display();

  // Sample the Model at regular steps
  VectorDouble hh = VH::sequence(0., 3., 3./50.);
  CovCalcMode mode(ECalcMember::LHS);
  mode.setAsVario(true);
  VH::display("\nModel sampled",modellmc.sample(hh,VectorDouble(),0,0,&mode));

  /////////////////////////////
  // Creating the Tapered Model
  CovLMCTapering covtape = CovLMCTapering(ETape::STORKEY, 4., ctxt.getSpace());
  // Build the Covariance list
  covtape.addCov(&cov1);
  covtape.addCov(&cov2);
  // Building the Model
  Model modeltape = Model(ctxt);
  modeltape.setCovList(&covtape);
  modeltape.display();

  // Sample the Tapered Model at regular steps
  VH::display("\nTapered Model",modeltape.sample(hh,VectorDouble(),0,0,&mode));

  /////////////////////////////
  // Creating the Convoluted Model
  CovLMCConvolution covconv = CovLMCConvolution(EConvType::EXPONENTIAL, EConvDir::X, 1., 10, ctxt.getSpace());
  // Build the Covariance list
  covconv.addCov(&cov1);
  covconv.addCov(&cov2);
  // Building the Model
  Model modelconv = Model(ctxt);
  modelconv.setCovList(&covconv);
  modelconv.display();
  // Sample the Tapered Model at regular steps
  VH::display("\nConvoluted Model", modelconv.sample(hh,VectorDouble(),0,0,&mode));

  /////////////////////////////////////////
  // Creating Covariance and Drift matrices
  Model* modelS = Model::createFromEnvironment(1, 3);
  modelS->addCovFromParam(ECov::CUBIC, 10., 12.);
  modelS->addCovFromParam(ECov::SPHERICAL, TEST, 23., TEST, {2., 3., 4.}, VectorDouble(),
                          {10., 20., 30.});
  modelS->addDrift(new DriftM());
  modelS->addDrift(new DriftM(VectorInt({1})));
  modelS->addDrift(new DriftM(VectorInt({2})));
  modelS->display();

  /////////////////////////////
  // Serialization of the Model
  modelS->dumpToNF("Complex");
  Model* modelSS = Model::createFromNF("Complex");
  modelSS->display();

  /////////////////////////////////////////////
  // Building the Covariance and Drift Matrices
  // (selection and heterotopy)
  MatrixSquareSymmetric covM;
  MatrixRectangular driftM;

  Model* modelM = Model::createFromEnvironment(2, 2);
  modelM->addCovFromParam(ECov::CUBIC, 10., TEST, 0., {}, {2., 1., 1., 4.});
  modelM->addDrift(new DriftM());                  // Universality Condition
  modelM->addDrift(new DriftM(VectorInt({1})));    // Drift: X
  modelM->addDrift(new DriftM(VectorInt({0,1})));  // Dirft: Y
  modelM->display();

  int nsample = workingDbc->getSampleNumber();
  // Adding a first variable (filled completely)
  VectorDouble rnd1 = VH::simulateGaussian(nsample);
  workingDbc->addColumns(rnd1, "Z1");
  // Adding a second variable (with one TEST values)
  VectorDouble rnd2 = VH::simulateGaussian(nsample);
  rnd2[1] = TEST;
  workingDbc->addColumns(rnd2, "Z2");
  VectorDouble verr1 = VectorDouble(nsample, 0.1);
  verr1[3] = TEST;
  workingDbc->addColumns(verr1, "V1");
  VectorDouble verr2 = VectorDouble(nsample, 0.25);
  workingDbc->addColumns(verr2, "V2");
  // Adding a Selection
  workingDbc->addColumns({1, 1, 1, 0, 1, 0}, "Sel");
  OptCst::define(ECst::NTCOL, -1);
  OptCst::define(ECst::NTROW, -1);
  DbStringFormat* dbfmt = DbStringFormat::createFromFlags(false, true, false, false, true);
  workingDbc->display(dbfmt);

  // Complete Matrices on the whole grid
  message("Covariance Matrix (complete)\n");
  MatrixSquareSymmetric covMS = modelM->evalCovMatrixSymmetric(workingDbc);
  covMS.display();

  message("Covariance Matrix Optimal (complete)\n");
  MatrixSquareSymmetric covMO = modelM->evalCovMatrixSymmetricOptim(workingDbc);
  covMO.display();

  message("Covariance Matrix Sparse (complete)\n");
  MatrixSparse* covMSS = modelM->evalCovMatrixSparse(workingDbc);
  covMSS->display();
  delete covMSS;

  message("Drift Matrix (complete)\n");
  driftM = modelM->evalDriftMatrix(workingDbc);
  driftM.display();

  // Adding the selection
  workingDbc->setLocator("Sel", ELoc::SEL);
  message("Covariance Matrix (with selection)\n");
  covM = modelM->evalCovMatrixSymmetric(workingDbc);
  covM.display();

  message("Drift Matrix (with selection)\n");
  driftM = modelM->evalDriftMatrix(workingDbc);
  driftM.display();

  // Adding the variables (heterotopic multivariate)
  workingDbc->setLocators({"Z*"}, ELoc::Z);
  message("Covariance Matrix (with selection & heterotopic multivariate)\n");
  covM = modelM->evalCovMatrixSymmetric(workingDbc);
  covM.display();

  message("Drift Matrix (with selection & heterotopic multivariate)\n");
  driftM = modelM->evalDriftMatrix(workingDbc);
  driftM.display();

  // Adding the variables (heterotopic multivariate & Verr)
  workingDbc->setLocators({"V*"}, ELoc::V);
  message("Covariance Matrix (with selection & heterotopic multivariate & verr)\n");
  covM = modelM->evalCovMatrixSymmetric(workingDbc);
  covM.display();

  message("Drift Matrix (with selection & heterotopic multivariate & verr)\n");
  driftM = modelM->evalDriftMatrix(workingDbc);
  driftM.display();

  // Selecting samples
  VectorInt nbgh = {0, 2, 3, 5};
  VH::display("Ranks of selected samples = ",nbgh);

  message("Covariance Matrix (selection & heterotopic multivariate & sampling)\n");
  covM = modelM->evalCovMatrixSymmetric(workingDbc, -1, nbgh);
  covM.display();

  message("Drift Matrix (selection & heterotopic multivariate & sampling)\n");
  driftM = modelM->evalDriftMatrix(workingDbc, -1, nbgh);
  driftM.display();

  // Testing Models on the Sphere

  defineDefaultSpace(ESpaceType::SN, 2);
  int ns = 20;
  int nincr = 30;
  VectorDouble incr = VH::sequence(0., GV_PI + EPSILON10, GV_PI / (nincr-1.));
  double mu = 1.0;
  double kappa = 2.0;

//  Model* modelSph = Model::createFromParam(ECov::LINEARSPH);
//  Model* modelSph = Model::createFromParam(ECov::GEOMETRIC, 0.9);
//  Model* modelSph = Model::createFromParam(ECov::POISSON, 1., 1., 10.);
//  Model* modelSph = Model::createFromParam(ECov::EXPONENTIAL, 5.0, 1., 0.,
//                                           VectorDouble(), VectorDouble(), VectorDouble(),
//                                           nullptr, true);
  Model *modelSph = Model::createFromParam(ECov::MATERN, 1./kappa, 1., mu,
                                           VectorDouble(), VectorDouble(),
                                           VectorDouble(), nullptr, false);
  VH::display("Spectrum", modelSph->getCova(0)->evalSpectrumOnSphere(ns));
  VH::display("Covariance", modelSph->getCova(0)->evalCovOnSphereVec(incr));

  return 0;
}
