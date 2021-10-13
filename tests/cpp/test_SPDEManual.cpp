#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "geoslib_e.h"

#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Model/NoStatFunctional.hpp"
#include "Mesh/AMesh.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/Law.hpp"
#include "Basic/FunctionalSpirale.hpp"

#include "MatrixC/MatrixCRectangular.hpp"

#include <algorithm>
#include <math.h>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "../../include/Mesh/MeshFactory.hpp"

#define __USE_MATH_DEFINES
#include <cmath>

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])

{
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("SPDEManual-");
  int seed = 10355;
  law_set_random_seed(seed);

  ///////////////////////
  // Creating the Db Grid
  auto nx={ 101,101 };
  Db workingDbc(nx);

  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);
  VectorDouble angle = spirale.getFunctionValues(&workingDbc);
  workingDbc.addFields(angle,"angle",ELoc::NOSTAT);

  //////////////////////
  //Creating the Mesh
  MeshETurbo mesh(workingDbc);

  ///////////////////////
  // Creating the Model
  Model model = Model(&workingDbc);
  CovAniso cova = CovAniso(ECov::BESSEL_K,model.getContext());
  cova.setRanges({10,45});
  model.addCova(&cova);

  /////////////////////////////////////////////////////
  // Creating the Precision Operator for simulation
  NoStatArray NoStat({"A"},&workingDbc);
  model.addNoStat(&NoStat);
  SPDE spde(model,workingDbc);

  ShiftOpCs S(&mesh, &model, &workingDbc);
  PrecisionOp Qsimu(&S, &cova, EPowerPT::MINUSHALF);

  ///////////////////////////
  // Simulation (Chebyshev)
  VectorDouble resultSimu;
  VectorDouble tab = ut_vector_simulate_gaussian(mesh.getNApices());

  resultSimu.resize(tab.size());
  Qsimu.eval(tab,resultSimu);
  workingDbc.addFields(resultSimu,"Simu",ELoc::Z);

  ///////////////////////////
  // Creating Data
  auto ndata = 1000;
  Db dat = Db(ndata, { 0., 0. }, { 100., 100. });

  ///////////////////////////
  // Simulating Data points
  ProjMatrix B(&dat, &mesh);
  VectorDouble datval(ndata);
  B.mesh2point(resultSimu, datval);
  dat.addFields(datval, "Simu", ELoc::Z);

  ///////////////////////////
  // Kriging
  double nug = 0.1;
  VectorDouble rhs(S.getSize());
  B.point2mesh(dat.getField("Simu"), rhs);
  for (auto &e : rhs)
  {
    e /= nug;
  }

  PrecisionOp Qkriging(&S, &cova, EPowerPT::ONE);
  PrecisionOpMultiConditional A;
  A.push_back(&Qkriging, &B);
  A.setVarianceData(0.01);

  VectorVectorDouble Rhs, resultvc;
  VectorDouble vc(S.getSize());

  resultvc.push_back(vc);
  Rhs.push_back(VectorDouble(rhs));

  A.evalInverse(Rhs, resultvc);
  workingDbc.addFields(resultvc[0], "Kriging");
  workingDbc.serialize("spde.ascii");
  return 0;
}
