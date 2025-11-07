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
#include "Basic/AStringable.hpp"
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

  // Identification
  message("On identifie le double\n");
  tab.identify();

  message("On identifie le int\n");
  itab.identify();

  // Basic manipulations
  tab.dump("Initial vector of Double");
  tab.add(tab);
  tab.dump("After adding this to itself");
  tab.addCst(12);
  tab.dump("After adding a constant to this");

  itab.dump("Initial Vector of Id");
  itab.add(itab);
  itab.dump("After adding this to itself");
  itab.addCst(12);
  itab.dump("After adding a constant to this");

  return 0;
}
