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

  mestitle(1, "Testing Interactive input");

  ianswer = questionId("Enter an Integer with no Default value");
  message("Value asentered interactivelyked = %d\n", ianswer);

  ianswer = questionId("Enter an Integer with Default value", 14);
  message("Value entered interactively = %d\n", ianswer);

  ranswer = questionDouble("Enter a Double with no Default value");
  message("Value entered interactively = %lf\n", ranswer);

  ranswer = questionDouble("Enter a Double with Default value", 14.);
  message("Value entered interactively = %lf\n", ranswer);

  mestitle(1, "Conversion from String to double/Id");

  auto dval = fromStrToDouble("  -12.345  ");
  message("Value read = %lf\n", dval);

  auto ival = fromStrToId("   12345 ");
  message("Value read = %d\n", ival);

  return 0;
}
