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

// This test is meant to demonstrate the Potential Model
// through estimation, cross-validation and simulations

#include "geoslib_d.h"
#include "geoslib_old_f.h"

#include "Enum/ECov.hpp"
#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCustom.hpp"
#include "Basic/Vector.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Potential-");

  // Global parameters
  int ndim = 2;
  ASpaceObject::defineDefaultSpace(ESpaceType::SPACE_RN, ndim);
  DbStringFormat dbfmt(FLAG_RESUME | FLAG_VARS | FLAG_ARRAY);

  // Generating the different constraining files

  // Iso-Potential file
  VectorDouble tabiso = { 7., 6., 1.,
                          5., 6., 1.,
                          6., 5., 1.,
                          3., 6., 2.,
                          7., 7., 2.,
                          8., 3., 2.,
                          8., 1., 3.,
                          7., 9., 3.,
                          10., 5., 3.,
                          3., 1., 3. };
  Db* dbiso = Db::createFromSamples(10, ELoadBy::SAMPLE, tabiso,
                                    {"x","y","iso"},
                                    {"x1","x2","layer"});
  dbiso->display(&dbfmt);

  // Gradient file
  VectorDouble tabgrd = { 1., 6., 1., 0.,
                          9., 2., -1., 1.,
                          7., 8., 0., -1 };
  Db* dbgrd = Db::createFromSamples(3, ELoadBy::SAMPLE, tabgrd,
                                    {"x","y","gx","gy"},
                                    {"x1","x2","g1","g2"});
  dbgrd->display(&dbfmt);

  // Tangent file
  VectorDouble tabtgt = { 3., 7., 1., 0.,
                          9., 7., 0.5, -0.5 };
  Db* dbtgt = Db::createFromSamples(2, ELoadBy::SAMPLE, tabtgt,
                                    {"x","y","tx","ty"},
                                    {"x1","x2","tangent1","tangent2"});
  dbtgt->display(&dbfmt);

  // Generate the output grid
  VectorInt nx = {101,101};
  VectorDouble dx = {0.1, 0.1};
  DbGrid* grid = DbGrid::create(nx, dx);
  grid->display();

  // Create the model
  Model* model = Model::createFromParam(ECov::CUBIC, 3.);

  // Create the Neighborhood (unique)
  NeighUnique* neighU = NeighUnique::create(ndim);

  // Launch the Potential estimation
  (void) potential_kriging(dbiso, dbgrd, dbtgt, grid, model, neighU,
                           0., 0., 0, 0, 0, 1);

//  (void) grid->dumpToNF("Grid.ascii");

  // ====================== Free pointers ==================================
  if (dbiso != nullptr) delete dbiso;
  if (grid != nullptr) delete grid;

  return (0);
}

