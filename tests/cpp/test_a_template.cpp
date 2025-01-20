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
#include "Db/DbGrid.hpp"
#include "Db/DbHelper.hpp"
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
  auto space1d = std::make_shared<const SpaceRN>(1);
  Model* mT = Model::createFromParam(ECov::EXPONENTIAL,
                                     scaleT,
                                     1.,
                                     1.,
                                     VectorDouble(),
                                     MatrixSquareSymmetric(),
                                     VectorDouble(),
                                     space1d,
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
                                    
 

  delete mT;
  delete mS;
  return(0);
}
