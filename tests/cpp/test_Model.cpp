/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/Law.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "Basic/File.hpp"
#include "Basic/VectorHelper.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Model-");
  int seed = 10355;
  int ndim = 2;
  int nvar = 1;
  law_set_random_seed(seed);

  ///////////////////////
  // Creating the Db
  auto nx={ 3,3 };
  DbGrid* workingDbc = DbGrid::create(nx);

  ///////////////////////
  // Creating the Model

  // Building the Covariance Context
  CovContext ctxt(nvar, ndim);
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
  // Building the Model
  Model modellmc = Model(ctxt);
  modellmc.setCovList(&covlmc);
  modellmc.display();

  ///////////////////////
  // Building the Covariance Matrix
  VectorDouble result = modellmc.covMatrixV(workingDbc, nullptr, 0, 0, 0, 1);

  // Checking that the matrix (VectorDouble) has been correctly filled by asking for statistics
  VH::displayStats("\nStatistics on Covariance Matrix",result);

  // Sample the Model at regular steps
  VectorDouble hh = VH::sequence(0., 3., 3./50.);
  VectorDouble vec1 = modellmc.sample(hh);
  VH::display("\nModel sampled",vec1);

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
  VectorDouble vec2 = modeltape.sample(hh);
  VH::display("\nTapered Model",vec2);

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
  VectorDouble vec3 = modelconv.sample(hh);
  VH::display("\nConvoluted Model", vec3);

  // Serialization of the Model
  Model* modelS = Model::createFromEnvironment(1, 3);
  modelS->addCovFromParam(ECov::CUBIC, 10., 12.);
  modelS->addCovFromParam(ECov::SPHERICAL, TEST, 23., TEST, {2., 3., 4.}, VectorDouble(),
                          {10., 20., 30.});
  modelS->addDrift(DriftFactory::createDriftFunc(EDrift::UC));
  modelS->addDrift(DriftFactory::createDriftFunc(EDrift::X));
  modelS->addDrift(DriftFactory::createDriftFunc(EDrift::F, CovContext(), 0));
  modelS->addDrift(DriftFactory::createDriftFunc(EDrift::F, CovContext(), 1));
  modelS->display();
  modelS->dumpToNF("Complex");

  Model* modelSS = Model::createFromNF("Complex");
  modelSS->display();

  delete workingDbc;
  return 0;
}

