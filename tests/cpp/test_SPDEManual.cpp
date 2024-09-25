/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Basic/Law.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "Basic/File.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Mesh/AMesh.hpp"
#include "Mesh/MeshETurbo.hpp"

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
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("test_SPDEManual-");
  int seed = 10355;
  law_set_random_seed(seed);

  ///////////////////////
  // Creating the Db Grid
  auto nx = { 101,101 };
  DbGrid* workingDbc = DbGrid::create(nx);

  //////////////////////
  //Creating the Mesh
  MeshETurbo mesh(workingDbc);

  ///////////////////////
  // Creating the Model
  Model* model = Model::createFromParam(ECov::MATERN, 1., 1., 1., {10., 45.});
  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);
  model->getCova(0)->makeAngleNoStatFunctional(&spirale);

  /////////////////////////////////////////////////
  // Creating the Precision Operator for simulation
  
  CovAniso* cova = model->getCova(0);
  ShiftOpCs S(&mesh, cova, workingDbc);
  PrecisionOp Qsimu(&S, cova);

  // /////////////////////////
  // // Simulation (Chebyshev)
  VectorDouble resultSimu = Qsimu.simulateOne();
  workingDbc->addColumns(resultSimu,"Simu",ELoc::Z);

  // ///////////////////////////
  // // Creating Data
  auto ndata = 1000;
  Db* dat = Db::createFromBox(ndata, workingDbc->getCoorMinimum(), workingDbc->getCoorMaximum(), 432432);

  // /////////////////////////
  // // Simulating Data points
  ProjMatrix B(dat, &mesh);
  VectorDouble datval(ndata);
  B.mesh2point(resultSimu, datval);
  dat->addColumns(datval, "Simu", ELoc::Z);

  // //////////
  // // Kriging
  double nug = 0.1;
  VectorDouble rhs(S.getSize());
  auto datv = dat->getColumn("Simu");
  B.point2mesh(datv, rhs);
  
  for (int i = 0; i < (int)rhs.size(); i++)
    rhs[i] /= nug;

  PrecisionOp Qkriging(&S, cova);
  PrecisionOpMultiConditional A;
  A.push_back(&Qkriging, &B);
  A.setVarianceData(0.01);

  std::vector<std::vector<double>> Rhs, resultvc;
  VectorDouble vc(S.getSize());

  resultvc.push_back(vc);
  Rhs.push_back(rhs);

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
