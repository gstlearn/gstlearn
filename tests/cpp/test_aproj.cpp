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
#include "Enum/ESPDECalcMode.hpp"

#include "Db/DbGrid.hpp"
#include "Basic/File.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "LinearOp/ProjMatrix.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 * It will be compiled but not run nor diff'ed.
 */
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("AAA_");

  // parameters of the final grid on [0,1]^2
  VectorInt nx({200,200});
  VectorDouble dx({1/199, 1/199});
  VectorDouble x0({0,0});

  // Grid meshing on [-0.25, 1.25]^2
  VectorInt nxm({100,100});
  VectorDouble dxm({1.5/99, 1.5/99});
  VectorDouble x0m({-0.25,-0.25});

  DbGrid* grid     = DbGrid::create(nx,dx,x0);
  DbGrid* gridExt  = DbGrid::create(nxm,dxm,x0m);
  MeshETurbo* mesh = new MeshETurbo(gridExt);

  MatrixSparse* ap = ProjMatrix(grid, mesh).getAproj();
  std::cout << ap->isSparse() << std::endl;

  return (0);
}
