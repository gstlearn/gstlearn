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

// This test is meant to demonstrate the Potential Model
// through estimation, cross-validation and simulations

#include "API/Potential.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovGradientAnalytic.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Enum/ECov.hpp"
#include "Enum/ESpaceType.hpp"
#include "Model/Model.hpp"
#include "Space/ASpaceObject.hpp"

using namespace gstlrn;

static ModelGeneric* st_duplicate_for_potential(const Model* model)
{
  auto* new_model = new ModelGeneric(*model->getContext());
  auto covp       = CovGradientAnalytic(*model->getCovAniso(0));
  new_model->setCov(&covp);

  DriftList* drifts = DriftFactory::createDriftListForGradients(model->getDriftList());
  new_model->setDriftList(drifts);

  return new_model;
}

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("test_Potential-");
  int mode     = 0;
  bool debug   = true;
  bool verbose = false;

  //============================================================//
  // Exemple in 1-D
  //============================================================//
  if (mode == 0 || mode == 1)
  {
    mestitle(0, "Example in 1-D");

    Id ndim = 1;
    defineDefaultSpace(ESpaceType::RN, ndim);

    // Iso-Potential file
    VectorDouble tabiso = {30., 1.,
                           80., 1.,
                           40., 2.,
                           50., 2.};
    Db* dbiso           = Db::createFromSamples(4, ELoadBy::SAMPLE, tabiso,
                                                {"x", "iso"},
                                                {"x1", "layer"});

    // Gradient file
    VectorDouble tabgrd = {0., 1.};
    Db* dbgrd           = Db::createFromSamples(1, ELoadBy::SAMPLE, tabgrd,
                                                {"x", "gx"},
                                                {"x1", "g1"});

    // Generate the output grid
    VectorInt nx    = {101};
    VectorDouble dx = {1};
    DbGrid* grid    = DbGrid::create(nx, dx);

    // Create the model
    double range            = scale2range(ECov::GAUSSIAN, 20.);
    Model* model            = Model::createFromParam(ECov::GAUSSIAN, range);
    ModelGeneric* new_model = st_duplicate_for_potential(model);

    if (debug) OptDbg::setReference(1);
    (void)krigingPotential(dbiso, dbgrd, nullptr, grid, new_model,
                           0., 0., true, true, false, true, 0, verbose);
    OptDbg::setReference(-1);

    // Visualize the results
    (void)grid->dumpToNF("Grid1D.NF");

    // Cross-validation
    (void)xvalidPotential(dbiso, dbgrd, nullptr, new_model,
                          0., 0., true, verbose);

    delete dbiso;
    delete dbgrd;
    delete grid;
    delete model;
    delete new_model;
  }

  //============================================================//
  // Exemple in 2-D
  //============================================================//
  if (mode == 0 || mode == 2)
  {
    mestitle(0, "Example in 2-D");

    Id ndim = 2;
    defineDefaultSpace(ESpaceType::RN, ndim);

    // Iso-Potential file
    VectorDouble tabiso = {7., 6., 1.,
                           5., 6., 1.,
                           6., 5., 1.,
                           3., 6., 2.,
                           7., 7., 2.,
                           8., 3., 2.,
                           8., 1., 3.,
                           7., 9., 3.,
                           10., 5., 3.,
                           3., 1., 3.};
    Db* dbiso           = Db::createFromSamples(10, ELoadBy::SAMPLE, tabiso,
                                                {"x", "y", "iso"},
                                                {"x1", "x2", "layer"});

    // Gradient file
    VectorDouble tabgrd = {1., 6., 1., 0.,
                           9., 2., -1., 1.,
                           7., 8., 0., -1};
    Db* dbgrd           = Db::createFromSamples(3, ELoadBy::SAMPLE, tabgrd,
                                                {"x", "y", "gx", "gy"},
                                                {"x1", "x2", "g1", "g2"});

    // Tangent file
    VectorDouble tabtgt = {3., 7., 1., 0.,
                           9., 7., 0.5, -0.5};
    Db* dbtgt           = Db::createFromSamples(2, ELoadBy::SAMPLE, tabtgt,
                                                {"x", "y", "tx", "ty"},
                                                {"x1", "x2", "tangent1", "tangent2"});

    // Generate the output grid
    VectorInt nx    = {101, 101};
    VectorDouble dx = {0.1, 0.1};
    DbGrid* grid    = DbGrid::create(nx, dx);

    // Create the model
    Model* model            = Model::createFromParam(ECov::CUBIC, 6.);
    ModelGeneric* new_model = st_duplicate_for_potential(model);
    (void)krigingPotential(dbiso, dbgrd, dbtgt, grid, new_model,
                           0., 0., true, false, false, false, 0, true);

    (void)grid->dumpToNF("Grid2D.NF");

    // Free pointers
    delete dbiso;
    delete dbgrd;
    delete dbtgt;
    delete grid;
    delete model;
    delete new_model;
  }

  return (0);
}
