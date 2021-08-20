#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "geoslib_e.h"
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
int main(int argc, char *argv[])

{
  auto pygst = std::string(std::getenv("PYGSTLEARN_DIR"));
  int seed = 10355;
  law_set_random_seed(seed);

  ///////////////////////
  // Creating the Db Grid
  auto nx={ 101,101 };
  Db workingDbc(nx);

  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);
  VectorDouble angle = spirale.getFunctionValues(&workingDbc);
  workingDbc.addFields(angle,"angle",LOC_NOSTAT);

  ///////////////////////
  // Creating the Model
  Model model = Model(&workingDbc);

  CovAniso cova = CovAniso(COV_BESSEL_K,model.getContext());
  cova.setRanges({20,20});
  model.addCova(&cova);

  model.addDrift({"1","f1"});

  NoStatArray NoStat({"A"},&workingDbc);
  model.addNoStat(&NoStat);

  model.display(1);

  ///////////////////////
  // Creating Data
  auto ndata = 500;
  Db dat = Db(ndata, { 0., 0. }, { 100., 100. });
  VectorDouble z = ut_vector_simulate_gaussian(ndata);
  VectorDouble drift = dat.getField("x.1");
  ut_vector_multiply_inplace(drift,0.1);
  ut_vector_add_inplace(z,drift);
  ut_vector_addval(z,10);
  dat.addFields(z,"variable",LOC_Z);
  dat.addFields(drift,"Drift",LOC_F);
  dat.display(1);

  ///////////////////////
  // Running SPDE
  SPDE spde(model,workingDbc,&dat,CALCUL_KRIGING);
  VectorDouble result = spde.computeCoeffs();
  ut_vector_display("Drift Coefficients:",result);
  return 0;
}

