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

// This test is mean to check the Block Kriging (essentially for printout)

#include "Basic/AStringFormat.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/File.hpp"
#include "Enum/ESpaceType.hpp"

#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "Model/Model.hpp"

#include "Estimation/CalcGlobal.hpp"

int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Global parameters
  int ndim = 2;
  law_set_random_seed(32131);
  AStringFormat format;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Parameters
  int nech        = 4;
  int nvar        = 1;
  bool flagSK     = false;

  // Generate the data base
  Db* data = Db::createFillRandom(nech, ndim, nvar, 0);
  data->setLocVariable(ELoc::Z, 1, 0, TEST);
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_ARRAY);
  data->display(dbfmt);

  // Generate the target file
  VectorInt nx = {5, 5};
  VectorDouble dx = {0.2, 0.2};
  DbGrid* target = DbGrid::create(nx, dx);

  // Create the Model
  double scale = 0.7;
  MatrixSquareSymmetric* sills =
    MatrixSquareSymmetric::createRandomDefinitePositive(nvar);
  Model* model = Model::createFromParam(ECov::EXPONENTIAL, scale, 0., 0., VectorDouble(),
                                        *sills, VectorDouble(), nullptr, false);
  if (flagSK)
  {
    VectorDouble means = VH::simulateGaussian(nvar);
    model->setMeans(means);
  }
  else
    model->setDriftIRF(1, 0);
  VectorDouble means(nvar, 0.);
  model->display();

  // Test on Global Kriging
  VectorInt ndisc = {3, 3};
  data->display();
  target->display();
  global_kriging(data, target, model, 0, true);

  // Free pointers

  delete sills;
  delete data;
  delete target;
  delete model;

  return (0);
}
