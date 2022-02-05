#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])

{
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("SPDEAPI-");
  int seed = 10355;
  int ndim = 2;
  int nvar = 1;
  law_set_random_seed(seed);

  ///////////////////////
  // Creating the Db
  auto nx={ 3,3 };
  DbGrid* workingDbc = DbGrid::create(nx);
  int nech = workingDbc->getSampleNumber(true);

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
  return(0);
  modellmc.display();

  ///////////////////////
  // Building the Covariance Matrix
  int nval = nech * nech * nvar * nvar;
  VectorDouble result(nval);
  modellmc.covMatrix(result, workingDbc, nullptr, 0, 0, 0, 1);

  // Checking that the matrix (VectorDouble) has been correctly filled by asking for statistics
  ut_vector_display_stats("\nStatistics on Covariance Matrix",result);

  // Sample the Model at regular steps
  VectorDouble vec1 = modellmc.sample(3., 50);
  ut_vector_display("\nModel sampled",vec1);

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
  VectorDouble vec2 = modeltape.sample(3., 50);
  ut_vector_display("\nTapered Model",vec2);

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
  VectorDouble vec3 = modelconv.sample(3., 50);
  ut_vector_display("\nConvoluted Model", vec3);

  delete workingDbc;
  return 0;
}

