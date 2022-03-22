#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
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
#include "Basic/File.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Mesh/MeshFactory.hpp"

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
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("SPDEManual-");
  int seed = 10355;
  law_set_random_seed(seed);

  ///////////////////////
  // Creating the Db Grid
  auto nx={ 101,101 };
  DbGrid* workingDbc = DbGrid::create(nx);

  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);
  VectorDouble angle = spirale.getFunctionValues(workingDbc);
  workingDbc->addColumns(angle,"angle",ELoc::NOSTAT);

  //////////////////////
  //Creating the Mesh
  MeshETurbo mesh(workingDbc);

  ///////////////////////
  // Creating the Model
  Model* model = Model::createFromDb(workingDbc);
  CovContext ctxt(model->getContext());
  CovLMC covs(ctxt.getSpace());
  CovAniso cova = CovAniso(ECov::BESSEL_K,ctxt);
  cova.setRanges({10,45});
  covs.addCov(&cova);
  model->setCovList(&covs);

  /////////////////////////////////////////////////////
  // Creating the Precision Operator for simulation
  NoStatArray NoStat({"A"},workingDbc);
  model->addNoStat(&NoStat);

  ShiftOpCs S(&mesh, model, workingDbc);
  PrecisionOp Qsimu(&S, &cova, EPowerPT::MINUSHALF);

  ///////////////////////////
  // Simulation (Chebyshev)
  VectorDouble resultSimu;
  VectorDouble tab = ut_vector_simulate_gaussian(mesh.getNApices());

  resultSimu.resize(tab.size());
  Qsimu.eval(tab,resultSimu);
  workingDbc->addColumns(resultSimu,"Simu",ELoc::Z);

  ///////////////////////////
  // Creating Data
  auto ndata = 1000;
  Db* dat = Db::createFromBox(ndata, {0.,0.}, {100.,100.});

  ///////////////////////////
  // Simulating Data points
  ProjMatrix B(dat, &mesh);
  VectorDouble datval(ndata);
  B.mesh2point(resultSimu, datval);
  dat->addColumns(datval, "Simu", ELoc::Z);

  ///////////////////////////
  // Kriging
  double nug = 0.1;
  VectorDouble rhs(S.getSize());
  B.point2mesh(dat->getColumn("Simu"), rhs);
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
  workingDbc->addColumns(resultvc[0], "Kriging");
  DbStringFormat dsf(FLAG_RESUME | FLAG_STATS);
  workingDbc->display(&dsf);
  (void) workingDbc->dumpToNF("spde.ascii");

  delete dat;
  delete workingDbc;
  delete model;
  return 0;
}
