#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Basic/FunctionalSpirale.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])

{
  ASerializable::setPrefixName("SPDEDrift-");
  int seed = 10355;
  law_set_random_seed(seed);

  ///////////////////////
  // Creating the Db Grid
  auto nx={ 101,101 };
  Db* workingDbc = Db::createFromGrid(nx);

  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);
  VectorDouble angle = spirale.getFunctionValues(workingDbc);
  workingDbc->addFields(angle,"angle",ELoc::NOSTAT);

  ///////////////////////
  // Creating the Model
  Model* model = Model::createFromDb(workingDbc);
  CovContext ctxt(model->getContext());
  CovLMC covs(ctxt.getSpace());
  CovAniso cova = CovAniso(ECov::BESSEL_K,ctxt);
  cova.setRanges({20,20});
  covs.addCov(&cova);
  model->setCovList(&covs);

  model->addDrift({"1","f1"});

  NoStatArray NoStat({"A"},workingDbc);
  model->addNoStat(&NoStat);

  model->display();

  ///////////////////////
  // Creating Data
  auto ndata = 500;
  Db* dat = Db::createFromBox(ndata, { 0., 0. }, { 100., 100. });
  VectorDouble z = ut_vector_simulate_gaussian(ndata);
  VectorDouble drift = dat->getField("x.1");
  ut_vector_multiply_inplace(drift,0.1);
  ut_vector_add_inplace(z,drift);
  ut_vector_addval(z,10);
  dat->addFields(z,"variable",ELoc::Z);
  dat->addFields(drift,"Drift",ELoc::F);
  dat->display();

  ///////////////////////
  // Running SPDE
  SPDE spde(model,workingDbc,dat,ESPDECalcMode::KRIGING);
  VectorDouble result = spde.computeCoeffs();
  ut_vector_display("Drift Coefficients:",result);

  delete dat;
  delete workingDbc;
  delete model;
  return 0;
}

