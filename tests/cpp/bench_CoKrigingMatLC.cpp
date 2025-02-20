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

// This test is mean to check the Cokriging for combined target variables

#include "Basic/AStringFormat.hpp"
#include "Basic/NamingConvention.hpp"
#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptCst.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Estimation/CalcKriging.hpp"

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
  OptCst::define(ECst::NTROW, -1);
  OptCst::define(ECst::NTCOL, -1);

  // Parameters
  bool verbose    = true;
  int nech        = 3;
  int nvar        = 2;
  bool flagSK     = false;
  bool flagUnique = false;

  // Generate the data base
  Db* data = Db::createFillRandom(nech, ndim, nvar, 0);
  data->setLocVariable(ELoc::Z, 1, 0, TEST);
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_ARRAY);
  data->display(dbfmt);

  // Generate the target file
  Db* target = Db::createFillRandom(1, ndim, 0, 0);

  // Create the Model
  int order = (flagSK) ? -1 : 0;
  Model* model = Model::createFillRandom(ndim, nvar, {ECov::EXPONENTIAL}, 1., order);
  model->display();
  

  // Neighborhood
  ANeigh* neigh;
  int nmaxi     = nech;
  double radius = 5.;
  if (flagUnique)
    neigh = NeighUnique::create();
  else
    neigh = NeighMoving::create(false, nmaxi, radius);

  // Define the verbose option
  if (verbose) OptDbg::setReference(1);

  // Test on Collocated CoKriging in Unique Neighborhood
  mestitle(0, "Without MATLC");
  kriging(data, target, model, neigh, EKrigOpt::POINT, true, true, false, VectorInt(),
          VectorInt());
  target->display(dbfmt);

  mestitle(0, "With MATLC");
  MatrixRectangular* matLC = MatrixRectangular::createFromVD({1, -3.}, 1, nvar);
  matLC->display();
  kriging(data, target, model, neigh, EKrigOpt::POINT, true, true, false, VectorInt(),
          VectorInt(), matLC);
  target->display(dbfmt);

  // Free pointers

  delete neigh;
  delete data;
  delete target;
  delete model;
  delete matLC;

  return (0);
}
