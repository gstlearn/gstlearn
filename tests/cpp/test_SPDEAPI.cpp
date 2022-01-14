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
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("SPDEAPI-");
  int seed = 10355;
  law_set_random_seed(seed);

  ///////////////////////
  // Creating the Db
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
  cova.setRanges({10,45});
  covs.addCov(&cova);
  model->setCovList(&covs);

  NoStatArray NoStat({"A"},workingDbc);
  model->addNoStat(&NoStat);
  model->display();

  ///////////////////////
  // Creating Data
  int ndata = 100;
  Db* dat = Db::createFromBox(ndata, { 0., 0. }, { 100., 100. });
  VectorDouble z = ut_vector_simulate_gaussian(ndata);
  dat->addFields(z,"variable",ELoc::Z);

  ///////////////////////
  // Running SPDE
  SPDE spde(model,workingDbc,dat,ESPDECalcMode::SIMUNONCOND);
//  SPDE spde(model,workingDbc,&dat,ESPDECalcMode::SIMUCOND);
  spde.compute();
  spde.query(workingDbc);
  workingDbc->serialize("spde_simunc.ascii");

  SPDE spde2(model,workingDbc,dat,ESPDECalcMode::SIMUCOND);
  spde2.compute();
  spde2.query(workingDbc);

  SPDE spde3(model,workingDbc,dat,ESPDECalcMode::KRIGING);
  spde3.compute();
  spde3.query(workingDbc);

  delete dat;
  delete workingDbc;
  delete model;
  return 0;
}

