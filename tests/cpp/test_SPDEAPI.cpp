#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "Basic/File.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Db/DbGrid.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"

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
  ASerializable::setPrefixName("SPDEAPI-");
  int seed = 10355;
  law_set_random_seed(seed);

  ///////////////////////
  // Creating the Db
  auto nx={ 101,101 };
  DbGrid* workingDbc = DbGrid::create(nx);

  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);
  VectorDouble angle = spirale.getFunctionValues(workingDbc);
  workingDbc->addColumns(angle,"angle",ELoc::NOSTAT);

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
  Db* dat = Db::createFromBox(ndata, {0.,0.}, {100.,100.}, 43246);
  VectorDouble z = VH::simulateGaussian(ndata);
  dat->addColumns(z,"variable",ELoc::Z);
  dat->display();

  ///////////////////////
  // Running SPDE
  SPDE spde(model,workingDbc,dat,ESPDECalcMode::SIMUNONCOND);
  //  SPDE spde(model,workingDbc,&dat,ESPDECalcMode::SIMUCOND);
  spde.compute();
  spde.query(workingDbc);
  (void) workingDbc->dumpToNF("spde_simunc.ascii");

  SPDE spde2(model,workingDbc,dat,ESPDECalcMode::SIMUCOND);
  spde2.compute();
  spde2.query(workingDbc);

  SPDE spde3(model,workingDbc,dat,ESPDECalcMode::KRIGING);
  spde3.compute();
  spde3.query(workingDbc);

  DbStringFormat dbfmt(FLAG_STATS,{"spde*"});
  workingDbc->display(&dbfmt);

  delete dat;
  delete workingDbc;
  delete model;
  return 0;
}

