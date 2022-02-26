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
#include "geoslib_d.h"
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"

/**
 * This file is meant to test the External Format Read /Write operations
 */
int main(int /*argc*/, char */*argv*/[])

{
  String filename;
  DbGrid* gridnew;

  // Global parameters
  int ndim = 2;
  int nvar = 1;

  setup_license("Demonstration");
  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Convert-");

  // Generate the output grid
  VectorInt nx = {30,25};
  DbGrid* grid = DbGrid::create(nx);
  DbStringFormat dbfmt(FLAG_STATS);

  // Create the Model
  CovContext ctxt(nvar);
  Model* model = Model::create(ctxt);
  CovLMC covs(ctxt.getSpace());
  CovAniso cova(ECov::SPHERICAL, 25., 0., 2., ctxt);
  covs.addCov(&cova);
  model->setCovList(&covs);

  // Perform a non-conditional simulation
  simtub(nullptr, grid, model);
  grid->display(&dbfmt);
  int icol = grid->getLastUID();
  VectorInt cols(1);
  cols[0] = icol;

  mestitle(1,"Writing VTK Grid");
  filename = ASerializable::buildFileName("VTK.grid");
  db_write_vtk(filename.c_str(), grid, cols);

  mestitle(1,"Writing Zycor Grid");
  filename = ASerializable::buildFileName("Zycor.grid");
  db_grid_write_zycor(filename.c_str(), grid, icol);

  mestitle(1,"Reading Zycor Grid");
  gridnew = db_grid_read_zycor(filename.c_str());
  gridnew->display(&dbfmt);

  mestitle(1,"Writing BMP Grid");
  filename = ASerializable::buildFileName("Bmp.grid");
  db_grid_write_bmp(filename.c_str(), grid, icol, 1, 1, 10);

  mestitle(1,"Reading BMP Grid");
  gridnew = db_grid_read_bmp(filename.c_str());
  gridnew->display(&dbfmt);

  mestitle(1,"Writing Irap Grid");
  filename = ASerializable::buildFileName("Irap.grid");
  db_grid_write_irap(filename.c_str(), grid, icol);

  mestitle(1,"Writing IfpEn Grid");
  filename = ASerializable::buildFileName("IfpEn.grid");
  db_grid_write_ifpen(filename.c_str(), grid, 1, &icol);

  mestitle(1,"Reading IfpEn Grid");
  gridnew = db_grid_read_ifpen(filename.c_str());
  gridnew->display(&dbfmt);

  mestitle(1,"Writing Eclipse Grid");
  filename = ASerializable::buildFileName("Eclipse.grid");
  db_grid_write_eclipse(filename.c_str(), grid, icol);

  mestitle(1,"Writing XYZ Grid");
  filename = ASerializable::buildFileName("XYZ.grid");
  db_grid_write_XYZ(filename.c_str(), grid, icol);

  mestitle(1,"Writing VTK Format");
  filename = ASerializable::buildFileName("VTK.file");
  db_write_vtk(filename.c_str(), grid, cols);

  mestitle(1,"Writing ArcGis Format");
  filename = ASerializable::buildFileName("ArcGis.grid");
  db_grid_write_arcgis(filename.c_str(), grid, icol);

  // Free the pointers
  delete grid;
  delete model;

  return (0);
}
