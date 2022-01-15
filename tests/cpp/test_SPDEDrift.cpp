#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Db/Db.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"


String getTestData(const String& filename)
{
  String exec_dir = ASerializable::getExecDirectory();
  // This path is compatible with CMake generation
  String filepath(exec_dir + "../../doc/data/Scotland/" + filename);

  return filepath;
}

int main(int /*argc*/, char */*argv*/[])
{
  bool verbose = true;
  String filepath = getTestData("temperatures.ascii");
  std::cout << filepath <<std::endl;
  Db* temperatures = Db::createFromNF(filepath,false,verbose);

  filepath = getTestData("grid.ascii");
  std::cout << filepath <<std::endl;
  Db* grid = Db::createFromNF(filepath,true,verbose);

  filepath = getTestData("model.ascii");
  std::cout << filepath <<std::endl;
  Model* model = Model::createFromNF(filepath,verbose);

  grid->display();
  temperatures->display();
  model->display();
  SPDE spde(model,grid,temperatures,ESPDECalcMode::KRIGING);

  spde.compute();
  spde.query(grid);

  filepath = getTestData("result.ascii");
  grid->dumpToNF(filepath,verbose);

  delete temperatures;
  delete grid;
  delete model;

  return 0;
}

