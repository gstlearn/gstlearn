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
/* This file is meant to demonstrate the managment of Rule object            */
/*                                                                            */
/******************************************************************************/
#include "Basic/File.hpp"
#include "Covariances/KernelMatern.hpp"
#include "LithoRule/Rule.hpp"
#include "Simulation/CalcSimuPGS.hpp"
#include "geoslib_f.h"

using namespace gstlrn;

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("test_PGS-");
  Id error = 0;
  Id ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  CovContext ctxt(1, 2, 1.);

  // Prepare the Discrete process with Discretized Option
  set_test_discrete(false);

  // Creating the Rule
  Rule* rule = Rule::createFromNames({"S", "T", "F1", "F2", "F3"});
  rule->setFaciesName(1, "MyFacies");
  rule->setFaciesValue(2, 6);

  rule->display();
  (void)rule->dumpToNF("truerule.NF");
  delete rule;

  message("Display Rule after Serialize/Deserialize\n");
  rule = Rule::createFromNF("truerule.NF");
  rule->display();
  delete rule;

  return static_cast<int>(error);
}
