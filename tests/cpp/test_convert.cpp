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
#include "Db/Db.hpp"
#include "Basic/File.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "OutputFormat/AOF.hpp"

/**
 * This file is meant to test the External Format Read /Write operations
 */
int main(int argc, char *argv[])

{
  String filename;
  DbGrid* gridnew = nullptr;

  // Global parameters
  int ndim = 2;
  int nvar = 1;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  defineDefaultSpace(ESpaceType::RN, ndim);
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Convert-");

  // Generate the output grid
  VectorInt nx = {30,25};
  DbGrid* grid = DbGrid::create(nx);
  DbStringFormat dbfmt(FLAG_STATS);

  // Create the Model
  CovContext ctxt(nvar);
  Model* model = Model::create(ctxt);
  CovAnisoList covs(ctxt);
  CovAniso cova(ECov::SPHERICAL, 25., 0., 2., ctxt);
  covs.addCov(&cova);
  model->setCovAnisoList(&covs);

  // Perform a non-conditional simulation
  simtub(nullptr, grid, model);
  grid->display(&dbfmt);
  int icol = grid->getLastUID();
  VectorInt cols(1);
  cols[0] = icol;

  mestitle(1,"Writing VTK Grid");
  filename = ASerializable::buildFileName(2, "VTK.grid");
  db_write_vtk(filename.c_str(), grid, cols);

  mestitle(1,"Writing Zycor Grid");
  filename = ASerializable::buildFileName(2, "Zycor.grid");
  db_grid_write_zycor(filename.c_str(), grid, icol);

  mestitle(1, "Reading Zycor Grid");
  delete gridnew;
  gridnew = db_grid_read_zycor(filename.c_str());
  gridnew->display(&dbfmt);

  mestitle(1,"Writing BMP Grid");
  filename = ASerializable::buildFileName(2, "Bmp.grid");
  db_grid_write_bmp(filename.c_str(), grid, icol, 1, 1, 10);

  mestitle(1, "Reading BMP Grid");
  delete gridnew;
  gridnew = db_grid_read_bmp(filename.c_str());
  gridnew->display(&dbfmt);

  mestitle(1,"Writing Irap Grid");
  filename = ASerializable::buildFileName(2, "Irap.grid");
  db_grid_write_irap(filename.c_str(), grid, icol);

  mestitle(1,"Writing IfpEn Grid");
  filename = ASerializable::buildFileName(2, "IfpEn.grid");
  db_grid_write_ifpen(filename.c_str(), grid, 1, &icol);

  mestitle(1, "Reading IfpEn Grid");
  delete gridnew;
  gridnew = db_grid_read_ifpen(filename.c_str());
  gridnew->display(&dbfmt);

  mestitle(1,"Writing Eclipse Grid");
  filename = ASerializable::buildFileName(2, "Eclipse.grid");
  db_grid_write_eclipse(filename.c_str(), grid, icol);

  mestitle(1,"Writing XYZ Grid");
  filename = ASerializable::buildFileName(2, "XYZ.grid");
  db_grid_write_XYZ(filename.c_str(), grid, icol);

  mestitle(1,"Writing VTK Format");
  filename = ASerializable::buildFileName(2, "VTK.file");
  db_write_vtk(filename.c_str(), grid, cols);

  mestitle(1,"Writing ArcGis Format");
  filename = ASerializable::buildFileName(2, "ArcGis.grid");
  db_grid_write_arcgis(filename.c_str(), grid, icol);

  // Free the pointers
  delete grid;
  delete model;
  delete gridnew;

  return (0);
}
