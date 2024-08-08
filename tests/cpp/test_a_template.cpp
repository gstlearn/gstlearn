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
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Space/ASpaceObject.hpp"
#include "Simulation/SimuSpectral.hpp"

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

  defineDefaultSpace(ESpaceType::SN);

  int ndim = 2;
  VectorInt nx = {360, 180};
  VectorDouble dx(2);
  for (int idim = 0; idim < ndim; idim++)
    dx[idim] = (double) nx[idim] / (nx[idim]-1) * GV_PI / 180.;
  DbGrid* grd = DbGrid::create(nx,dx,{0,0});
  (void) grd->setName("x1", "phi");
  (void) grd->setName("x2", "theta");
  grd->display();

  int nd = 100;
  int ns = 100; // 10000;
  int seed = 132674;

  String model_type = "POISSON";
  Model* modelSph = Model::createFromParam(ECov::POISSON, 1., 1., 10.);

  SimuSpectral sim(modelSph);
  sim.simulateOnSphere(ns, nd, seed, false);
  sim.computeOnSphere(grd, false);

  grd->display();
  return(0);
}
