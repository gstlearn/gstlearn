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
/**
 * This file is meant to parametrized the ModelGeneric in terms of ParamInfo
 * and to fit the values of these parameters according to the Maximum LogLikelihood
 * method and using the Vecchia approximation.
 */
#include "Covariances/CovAniso.hpp"
#include "Enum/ESpaceType.hpp"
#include "Model/Model.hpp"
#include "Space/ASpaceObject.hpp"
#include "geoslib_define.h"
using namespace gstlrn;

int main(int argc, char* argv[])
{
  DECLARE_UNUSED(argc);
  DECLARE_UNUSED(argv);

  defineDefaultSpace(ESpaceType::RN, 3);

  double range = 10.;
  Model* model = Model::createFromParam(ECov::LINEAR, range);
  model->display();

  model->getCovAniso(0)->getCorAnisoModify()->setRanges({30, 20, 10});
  model->display();

  exit(0);
}
