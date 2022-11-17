#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This test is meant to check the Automatic Model Fitting facility
 ** for various Space Dimensions and various Calculation criteria.
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("AutoFit-");
  int seed = 10355;
  law_set_random_seed(seed);

  ////
  //// Basic 2-D example
  ////

  // Defining the Space Dimension
  int ndim = 2;
  ASpaceObject::defineDefaultSpace(ESpaceType::SPACE_RN, ndim);
  mestitle(0,"Testing Model Fitting in 2-D");

  // Defining a Model for simulating a data set
  Model* model = Model::createFromParam(ECov::CUBIC, 20.);

  // Defining a Data Base
  Db* db = Db::createFromBox(100, {0.,0.}, {100., 100.});

  // Simulate a Gaussian Random Function on the Data Base
  (void) simtub(nullptr, db, model);

  // Calculate the experimental variogram
  VarioParam* varioparam = VarioParam::createOmniDirection(10);
  Vario* vario = Vario::create(varioparam, db);
  vario->compute();
  vario->display();
  vario->dumpToNF("Vario2D");

  // Fitting an omni-directional model
  Model* model_fit = Model::createFromEnvironment(1, ndim);
  model_fit->fit(vario, {ECov::GAUSSIAN, ECov::LINEAR});
  model_fit->display();
  model_fit->dumpToNF("Model2D");

  delete model;
  delete db;
  delete varioparam;
  delete vario;
  delete model_fit;

  ////
  //// Basic 4-D example
  ////

  // Defining the Space Dimension
  ndim = 4;
  ASpaceObject::defineDefaultSpace(ESpaceType::SPACE_RN, ndim);
  mestitle(0,"Testing Model Fitting in 4-D");

  // Defining a Data Base
  db = Db::createFromBox(100, {0.,0.,0.,0.}, {100., 100., 100., 100.});
  VectorDouble tab = VH::simulateGaussian(db->getActiveSampleNumber(), 0., 1.);
  db->addColumns(tab, "Var", ELoc::Z);

  // Calculate the experimental variogram
  varioparam = VarioParam::createOmniDirection(20);
  vario = Vario::create(varioparam, db);
  vario->compute();
  vario->display();
  vario->dumpToNF("Vario4D");

  // Fitting an omni-directional model
  model_fit = Model::createFromEnvironment(1, ndim);
  model_fit->fit(vario, {ECov::GAUSSIAN, ECov::LINEAR});
  model_fit->display();
  model_fit->dumpToNF("Model4D");

  delete db;
  delete varioparam;
  delete vario;
  delete model_fit;

  return 0;
}

