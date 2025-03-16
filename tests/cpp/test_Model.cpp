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
#include "Covariances/CovAnisoList.hpp"
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
#include "Matrix/MatrixSquareSymmetric.hpp"

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
  CovAnisoList covlmc = CovAnisoList(ctxt);
  // Build the Elementary Covariances
  CovAniso cov1 = CovAniso(ECov::CUBIC,ctxt);
  cov1.setRanges({1.2,2.1});
  cov1.setSill(1.5);
  covlmc.addCov(&cov1);
  CovAniso cov2 = CovAniso(ECov::NUGGET,ctxt);
  cov2.setSill(0.5);
  covlmc.addCov(&cov2);
  // Assembling the Model
  modellmc.setCovAnisoList(&covlmc);
  modellmc.display();

  // Building the Covariance Matrix
  MatrixSquareSymmetric result = modellmc.evalCovMatSym(workingDbc);
  result.display();

  // Sample the Model at regular steps
  VectorDouble hh = VH::sequence(0., 3., 3./50.);
  CovCalcMode mode(ECalcMember::LHS);
  mode.setAsVario(true);
  VH::dump("\nModel sampled",modellmc.sample(hh,VectorDouble(),0,0,&mode));

  /////////////////////////////
  // Creating the Tapered Model
  CovLMCTapering covtape = CovLMCTapering(ETape::STORKEY, 4., ctxt);
  // Build the Covariance list
  covtape.addCov(&cov1);
  covtape.addCov(&cov2);
  // Building the Model
  Model modeltape = Model(ctxt);
  modeltape.setCovAnisoList(&covtape);
  modeltape.display();

  // Sample the Tapered Model at regular steps
  VH::dump("\nTapered Model",modeltape.sample(hh,VectorDouble(),0,0,&mode));

  /////////////////////////////
  // Creating the Convoluted Model
  CovLMCConvolution covconv = CovLMCConvolution(EConvType::EXPONENTIAL, EConvDir::X, 1., 10, ctxt);
  // Build the Covariance list
  covconv.addCov(&cov1);
  covconv.addCov(&cov2);
  // Building the Model
  Model modelconv = Model(ctxt);
  modelconv.setCovAnisoList(&covconv);
  modelconv.display();
  // Sample the Tapered Model at regular steps
  VH::dump("\nConvoluted Model", modelconv.sample(hh,VectorDouble(),0,0,&mode));

  /////////////////////////////////////////
  // Creating Covariance and Drift matrices
  Model* modelS = Model::createFromEnvironment(1, 3);
  modelS->addCovFromParam(ECov::CUBIC, 10., 12.);
  modelS->addCovFromParam(ECov::SPHERICAL, TEST, 23., TEST, {2., 3., 4.},
                          MatrixSquareSymmetric(), {10., 20., 30.});
  DriftM FF;
  modelS->addDrift(&FF);
  FF = DriftM(VectorInt({1}));
  modelS->addDrift(&FF);
  FF = DriftM(VectorInt({2}));
  modelS->addDrift(&FF);
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
  MatrixSquareSymmetric* sills =
    MatrixSquareSymmetric::createFromVD({2., 1., 1., 4.});
  modelM->addCovFromParam(ECov::CUBIC, 10., TEST, 0., VectorDouble(), *sills);
  delete sills;
  FF = DriftM(); // Universality Condition
  modelM->addDrift(&FF);
  FF = DriftM(VectorInt({1})); // Drift: X
  modelM->addDrift(&FF); 
  FF = DriftM(VectorInt({0,1})); // Drift: Y
  modelM->addDrift(&FF);
  modelM->display();

  int nsample = workingDbc->getNSample();
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
  MatrixSquareSymmetric covMS = modelM->evalCovMatSym(workingDbc);
  covMS.display();

  message("Covariance Matrix Optimal (complete)\n");
  MatrixSquareSymmetric covMO = modelM->evalCovMatSym(workingDbc);
  covMO.display();

  message("Covariance Matrix Sparse (complete)\n");
  MatrixSparse* covMSS = modelM->evalCovMatSparse(workingDbc);
  covMSS->display();
  delete covMSS;

  message("Drift Matrix (complete)\n");
  driftM = modelM->evalDriftMat(workingDbc);
  driftM.display();

  // Adding the selection
  workingDbc->setLocator("Sel", ELoc::SEL, 0);
  message("Covariance Matrix (with selection)\n");
  covM = modelM->evalCovMatSym(workingDbc);
  covM.display();

  message("Drift Matrix (with selection)\n");
  driftM = modelM->evalDriftMat(workingDbc);
  driftM.display();

  // Adding the variables (heterotopic multivariate)
  workingDbc->setLocators({"Z*"}, ELoc::Z, 0);
  message("Covariance Matrix (with selection & heterotopic multivariate)\n");
  covM = modelM->evalCovMatSym(workingDbc);
  covM.display();

  message("Drift Matrix (with selection & heterotopic multivariate)\n");
  driftM = modelM->evalDriftMat(workingDbc);
  driftM.display();

  // Adding the variables (heterotopic multivariate & Verr)
  workingDbc->setLocators({"V*"}, ELoc::V, 0);
  message("Covariance Matrix (with selection & heterotopic multivariate & verr)\n");
  covM = modelM->evalCovMatSym(workingDbc);
  covM.display();

  message("Drift Matrix (with selection & heterotopic multivariate & verr)\n");
  driftM = modelM->evalDriftMat(workingDbc);
  driftM.display();

  // Selecting samples
  VectorInt nbgh = {0, 2, 3, 5};
  VH::dump("Ranks of selected samples = ",nbgh);

  message("Covariance Matrix (selection & heterotopic multivariate & sampling)\n");
  covM = modelM->evalCovMatSym(workingDbc, nbgh, -1);
  covM.display();

  message("Drift Matrix (selection & heterotopic multivariate & sampling)\n");
  driftM = modelM->evalDriftMat(workingDbc, nbgh);
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
//                                           VectorDouble(), MatrixSquareSymmetric(), VectorDouble(),
//                                           nullptr, true);
  Model *modelSph = Model::createFromParam(ECov::MATERN, 1./kappa, 1., mu,
                                           VectorDouble(), MatrixSquareSymmetric(),
                                           VectorDouble(), nullptr, false);
  VH::dump("Spectrum", modelSph->getCovAniso(0)->evalSpectrumOnSphere(ns));
  VH::dump("Covariance", modelSph->getCovAniso(0)->evalCovOnSphereVec(incr));


  delete workingDbc;
  delete modelM;
  delete modelSS;
  delete modelS;
  delete modelSph;
  delete dbfmt;

  return 0;
}
