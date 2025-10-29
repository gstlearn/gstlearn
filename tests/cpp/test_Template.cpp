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

  // Use the traditional function cumul() in its real version (polymorphism)
  std::cout << "cumul double = " << VH::cumul(tab) << std::endl;
  // Use the traditional function cumul() in its integer version (polymorphism)
  std::cout << "cumul entier = " << VH::cumul(itab) << std::endl;
  // Use the template function cumul2() in its real version
  std::cout << "cumul neutre (tab) = " << VH::cumul2(tab) << std::endl;
  // Use the template function cumul2() in its integer version
  std::cout << "cumul neutre (itab) = " << VH::cumul2(itab) << std::endl;

  return 0;
}
