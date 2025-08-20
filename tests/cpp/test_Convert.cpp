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
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "OutputFormat/AOF.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

using namespace gstlrn;

/**
 * This file is meant to test the External Format Read /Write operations
 */
int main(int argc, char* argv[])

{
  String filename;
  DbGrid* gridnew = nullptr;

  // Global parameters
  Id ndim = 2;
  Id mode = 0;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  defineDefaultSpace(ESpaceType::RN, ndim);
  ASerializable::setPrefixName("test_Convert-");

  // Generate the output grid
  VectorInt nx = {30, 25};
  DbGrid* grid = DbGrid::create(nx);
  DbStringFormat dbfmt(FLAG_STATS);

  // Create the Model
  Model* model = Model::createFromParam(ECov::SPHERICAL, 25., 0., 2.);

  // Perform a non-conditional simulation
  simtub(nullptr, grid, model);
  grid->display(&dbfmt);
  auto icol = grid->getLastUID();
  VectorInt cols(1);
  cols[0] = icol;

  // VTK
  if (mode == 0 || mode == 1)
  {
    mestitle(1, "Writing VTK Grid");
    filename = ASerializable::buildFileName(2, "VTK.grid");
    db_write_vtk(filename.c_str(), grid, cols);
  }

  // ZYCOR
  if (mode == 0 || mode == 2)
  {
    mestitle(1, "Writing Zycor Grid");
    filename = ASerializable::buildFileName(2, "Zycor.grid");
    db_grid_write_zycor(filename.c_str(), grid, icol);

    mestitle(1, "Reading Zycor Grid");
    delete gridnew;
    gridnew = db_grid_read_zycor(filename.c_str());
    gridnew->display(&dbfmt);
  }

  // BMP
  if (mode == 0 || mode == 3)
  {
    mestitle(1, "Writing BMP Grid");
    filename = ASerializable::buildFileName(2, "Bmp.grid");
    db_grid_write_bmp(filename.c_str(), grid, icol, 1, 1, 10);

    mestitle(1, "Reading BMP Grid");
    delete gridnew;
    gridnew = db_grid_read_bmp(filename.c_str());
    gridnew->display(&dbfmt);
  }

  // IRAP
  if (mode == 0 || mode == 4)
  {
    mestitle(1, "Writing IRAP Grid");
    filename = ASerializable::buildFileName(2, "Irap.grid");
    db_grid_write_irap(filename.c_str(), grid, icol);
  }

  // IFPEN
  if (mode == 0 || mode == 5)
  {
    mestitle(1, "Writing IfpEn Grid");
    filename = ASerializable::buildFileName(2, "IfpEn.grid");
    db_grid_write_ifpen(filename.c_str(), grid, 1, &icol);

    mestitle(1, "Reading IfpEn Grid");
    delete gridnew;
    gridnew = db_grid_read_ifpen(filename.c_str());
    gridnew->display(&dbfmt);
  }

  // Eclipse
  if (mode == 0 || mode == 6)
  {
    mestitle(1, "Writing Eclipse Grid");
    filename = ASerializable::buildFileName(2, "Eclipse.grid");
    db_grid_write_eclipse(filename.c_str(), grid, icol);
  }

  // XYZ
  if (mode == 0 || mode == 7)
  {
    mestitle(1, "Writing XYZ Grid");
    filename = ASerializable::buildFileName(2, "XYZ.grid");
    db_grid_write_XYZ(filename.c_str(), grid, icol);
  }

  // ArcGis
  if (mode == 0 || mode == 8)
  {
    mestitle(1, "Writing ArcGis Format");
    filename = ASerializable::buildFileName(2, "ArcGis.grid");
    db_grid_write_arcgis(filename.c_str(), grid, icol);
  }

  // Free the pointers
  delete grid;
  delete model;
  delete gridnew;

  return (0);
}
