/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Db/Db.hpp"
#include "Basic/File.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Variogram/DirParam.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Geometry/Geometry.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 */
int main(int /*argc*/, char */*argv*/[])

{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  int ndim = 1;
  ASpaceObject::defineDefaultSpace(ESpaceType::SPACE_RN, ndim);

  int num_steps=16;
  Db* mydb_real = new Db();
  VectorDouble xvec = ut_vector_sequence(1., num_steps, 1.);
  mydb_real->addColumns(xvec, "x", ELoc::X);
  VectorDouble zvec = {0.77745871, 0.77786762, 0.80903279, 0.80464491, 0.73316259, 0.73454369, 0.83387421, 0.74832878, 0.65179765, 0.65792214, 0.6767367, 0.68796146, 0.69657135, 0.62711037, 0.55293927, 0.50427098};
  mydb_real->addColumns(zvec, "z", ELoc::Z);
  mydb_real->display();

  DirParam* mydir = new DirParam(1,20,1);
  VarioParam* myVarioParamOmni = new VarioParam();
  myVarioParamOmni->addDir(*mydir);
  Vario* myVarioOmni = Vario::create(myVarioParamOmni,mydb_real);
  myVarioOmni->compute(ECalcVario::VARIOGRAM);
  myVarioOmni->display();
  Model* mymodel = Model::createFromDb(mydb_real);

  (void) mymodel->fit(myVarioOmni,{ECov::EXPONENTIAL});
  mymodel->display();

  return (0);
}
