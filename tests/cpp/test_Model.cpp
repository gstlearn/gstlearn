#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

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
  auto nx={ 3,3 };
  Db workingDbc(nx);


  ///////////////////////
  // Creating the Model
  Model model = Model(&workingDbc);
  CovAniso cova = CovAniso(ECov::CUBIC,model.getContext());
  cova.setRanges({10,45});
  model.addCova(&cova);


  model.display(1);

  VectorDouble result(81);

  model_covmat(&model,&workingDbc,nullptr,0,0,0,1,result.data());
  //model.covMatrix(&workingDbc,nullptr,0,0,0,1,result);
  return 0;
}

