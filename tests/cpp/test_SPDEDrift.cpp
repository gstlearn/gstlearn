#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Db/Db.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Basic/ASerializable.hpp"

int main(int /*argc*/, char */*argv*/[])
{
  bool verbose = true;
  String filename;

  filename = ASerializable::getTestData("Scotland","temperatures.ascii");
  std::cout << filename << std::endl;
  Db* temperatures = Db::createFromNF(filename,false,verbose);

  filename = ASerializable::getTestData("Scotland","grid.ascii");
  std::cout << filename <<std::endl;
  Db* grid = Db::createFromNF(filename,true,verbose);

  filename = ASerializable::getTestData("Scotland","model.ascii");
  std::cout << filename <<std::endl;
  Model* model = Model::createFromNF(filename,verbose);

  grid->display();
  temperatures->display();
  model->display();

  SPDE spde(model,grid,temperatures,ESPDECalcMode::KRIGING);

  spde.compute();
  spde.query(grid);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Drift-");
  grid->dumpToNF("result.ascii",verbose);

  delete temperatures;
  delete grid;
  delete model;

  return 0;
}

