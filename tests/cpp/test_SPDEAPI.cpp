#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "geoslib_e.h"

#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "Basic/Law.hpp"

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
  // Création de la db //

  auto nx={ 101,101 };
  Db workingDbc(nx);

  //////////////////////
  //Creating the Mesh
  MeshETurbo mesh(workingDbc);

  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);
  VectorDouble angle = spirale.getFunctionValues(&workingDbc);
  workingDbc.addFields(angle,"angle",LOC_NOSTAT);

  ///////////////////////
  // Creating the Model
  Model model = Model(&workingDbc);
  CovAniso cova = CovAniso(COV_BESSEL_K,model.getContext());
  cova.setRanges({10,45});
  model.addCova(&cova);


  NoStatArray NoStat({"A"},&workingDbc);
  model.addNoStat(&NoStat);

  // Creating Data
  auto ndata = 1000;
  Db dat = Db(ndata, { 0., 0. }, { 100., 100. });

  SPDE spde(model,workingDbc,&dat,CALCUL_KRIGING);
  spde.computeKriging();
  spde.query(&workingDbc);
  workingDbc.serialize(pygst + "spde.ascii");
  return 0;
}
