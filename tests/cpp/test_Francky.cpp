#include "geoslib_f.h"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Model/NoStatFunctional.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighUnique.hpp"

#include <algorithm>
#include <math.h>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#define __USE_MATH_DEFINES
#include <cmath>

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])

{
  int seed = 10355;
  law_set_random_seed(seed);

  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Francky-");
  // Creating the 2-D Db
  auto nx = { 101, 101 };
  DbGrid* workingDbc = DbGrid::create(nx);

  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);

  // Creating the Non-stationary Model
  Model model = Model(workingDbc);
  CovContext ctxt = model.getContext();
  CovLMC covs(ctxt.getSpace());
  CovAniso cova = CovAniso(ECov::BESSEL_K,ctxt);
  cova.setRanges({10,45});
  covs.addCov(&cova);
  model.setCovList(&covs);

  NoStatFunctional NoStat(&spirale);
  model.addNoStat(&NoStat);

  // Creating the 2-D Data Db with a Normal Variable
  auto ndata = 100;
  Db* dat = Db::createFromBox(ndata, {0.,0.}, {100.,100.});
  VectorDouble Z = ut_vector_simulate_gaussian(ndata);
  dat->addColumns(Z, "Z",ELoc::Z);

  // Creating the Neighborhood (Unique)
  NeighUnique* neighU = NeighUnique::create(2,false);

  // Testing Kriging
  kriging(dat,workingDbc,&model,neighU);
  (void) workingDbc->dumpToNF("franckyFunctional.ascii");

  DbStringFormat dbfmt(FLAG_STATS,{"Kriging*"});
  workingDbc->display(&dbfmt);

  message("Test performed successfully\n");

  delete dat;
  delete workingDbc;
  delete neighU;
  return 0;
}
