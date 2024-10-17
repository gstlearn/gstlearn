/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/*                                                                            */
/* This file is meant to demonstrate the PCA feature                          */
/*                                                                            */
/******************************************************************************/
#include "Enum/ECov.hpp"

#include "Enum/EDbg.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"
#include "Basic/OptDbg.hpp"
#include "geoslib_define.h"

/****************************************************************************/
/*!
** Main Program for testing MAF
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  DECLARE_UNUSED(argc, argv)
  int nlag = 10;
  double lag = 0.025;
  VarioParam* varioparam = VarioParam::createOmniDirection(nlag, lag);

  String filename;
  filename = ASerializable::getTestData("MAF", "data_for_MAF.dat");
  Db* data = Db::createFromNF(filename);
  data->display();

  Vario* vario_MAF = Vario::computeFromDb(*varioparam, data);
  vario_MAF->display();

  Model* model_MAF = Model::create();
  Constraints ctr;
  Option_VarioFit ovf;
  Option_AutoFit oaf;
  oaf.setVerbose(true);
  auto types = ECov::fromKeys({"NUGGET", "EXPONENTIAL", "SPHERICAL"});
  bool verbose = true;
  OptDbg::define(EDbg::CONVERGE);
  int err = model_MAF->fit(vario_MAF, types, ctr, ovf, oaf, verbose);
  model_MAF->display();

  delete varioparam;
  delete data;
  delete vario_MAF;
  delete model_MAF;

  return (err);
}
