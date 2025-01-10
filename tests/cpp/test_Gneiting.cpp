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
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CorGneiting.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceRN.hpp"
/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 * It will be compiled but not run nor diff'ed.
 */
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  VectorDouble scales = {2.,3.};
  double scaleT = 5.3;
  VectorDouble coords1 = {12.,3.,1.};
  VectorDouble coords2 = {4.,5.,2.};
  auto space1d = SpaceRN(1);
  double sep = 1.;
  Model* mT = Model::createFromParam(ECov::EXPONENTIAL,
                                     scaleT,
                                     1.,
                                     1.,
                                     VectorDouble(),
                                     MatrixSquareSymmetric(),
                                     VectorDouble(),
                                     &space1d,
                                     false);
                                    
  Model* mS = Model::createFromParam(ECov::EXPONENTIAL,
                                     1.,
                                     1.,
                                     1.,
                                     scales,
                                     MatrixSquareSymmetric(),
                                     VectorDouble(),
                                     nullptr,
                                     false);
                                    
  CovAniso* covT = mT->getCova(0);
  CovAniso* covS = mS->getCova(0);
  CorGneiting gneiting = CorGneiting(covS,covT,sep);
  SpacePoint p1(gneiting.getSpace());
  SpacePoint p2(gneiting.getSpace());
  p1.setCoords(coords1);
  p2.setCoords(coords2);
  double cres = gneiting.eval(p1,p2);
  std::cout << "Value of Gneiting " << cres <<std::endl;
  delete mT;
  delete mS;
  return(0);
}
