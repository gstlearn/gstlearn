#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Law.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "geoslib_e.h"

#include "geoslib_f.h"

#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Model/NoStatFunctional.hpp"
#include "Basic/FunctionalSpirale.hpp"

#include "MatrixC/MatrixCRectangular.hpp"

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
int main(int argc, char *argv[])

{
  int seed = 10355;
  law_set_random_seed(seed);
  auto pygst = std::string(std::getenv("PYGSTLEARN_DIR"));

  // Creating the 2-D Db
  auto nx = { 101, 101 };
  Db workingDbc(nx);

  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);

  // Creating the Non-stationary Model
  Model model = Model(&workingDbc);
  CovAniso cova = CovAniso(COV_BESSEL_K,model.getContext());
  cova.setRanges({10,45});
  model.addCova(&cova);
  NoStatFunctional NoStat(&spirale);
  model.addNoStat(&NoStat);

  // Creating the 2-D Data Db with a Normal Variable
  auto ndata = 100;
  Db dat = Db(ndata, {0.,0.}, {100.,100.});
  VectorDouble Z = ut_vector_simulate_gaussian(ndata);
  dat.addFields(Z, "Z",LOC_Z);

  // Creating the Neighborhood (Unique)
  Neigh neigh(2);

  // Testing Kriging
  kriging(&dat,&workingDbc,&model,&neigh);
  workingDbc.serialize(pygst+ "franckyFunctional.ascii");

  return 0;
}
