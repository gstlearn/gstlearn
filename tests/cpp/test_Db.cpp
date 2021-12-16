#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of the Db
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])

{
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("TestDb-");
  int seed = 10355;
  law_set_random_seed(seed);

  // Creating the Grid Rotated Db
  Db grid({6,4}, {1.,2.}, {10.,20.}, {10.,0.});
  grid.display();

  // Creating the Model
  Model model = Model(&grid);
  CovAniso cova = CovAniso(ECov::CUBIC,model.getContext());
  cova.setRanges({10,45});
  cova.setAnisoAngles({30.,0.});
  model.addCova(&cova);
  model.display(1);

  // Creating the MeshTurbo which contains the Db
  MeshETurbo mesh;
  mesh.initFromCova(cova,grid,10,2,true,1);

  return 0;
}

