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
#include "Basic/VectorHelper.hpp"

using namespace gstlrn;
/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to test the Template
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Establish a sequence of real values
  VectorDouble tab = VH::sequenceVD(0, 9., 1.);
  tab[3]           = TEST;

  // Establish a sequence of integer values
  VectorInt itab = VH::sequence(10);
  itab[3]        = ITEST;

  message("On identifie le double\n");
  tab.identify();

  message("On identifie le int\n");
  itab.identify();

  return 0;
}
