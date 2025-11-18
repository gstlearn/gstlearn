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
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"

/**
 * This test is meant to check the interactive questioning
 */
using namespace gstlrn;
/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main()
{
  Id ianswer;
  double ranswer;

  // Testing numerical input

  message("Testing Interactive input\n");

  ianswer = askInteractive<Id>("Enter an Integer with no Default value", getNA<Id>());
  message("Value read = %d\n", ianswer);

  ianswer = askInteractive<Id>("Enter an Integer with Default value", 14);
  message("Value read = %d\n", ianswer);

  ranswer = askInteractive<double>("Enter a Double with no Default value", getNA<double>());
  message("Value read = %lf\n", ranswer);

  ranswer = askInteractive<double>("Enter a Double with Default value", 14.);
  message("Value read = %lf\n", ranswer);

  return 0;
}
