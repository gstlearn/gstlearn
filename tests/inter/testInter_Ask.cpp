/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Basic/String.hpp"
#include "Basic/AStringable.hpp"
#include <string.h>

/**
 * This test is meant to check the interactive questioning
 * @return
 */
int main()
{
  int ianswer;
  double ranswer;

  // Testing numerical input

  message("Testing Interactive input\n");

  ianswer = askInt("Enter an Integer with no Default value");
  message("Value read = %d\n",ianswer);

  ianswer = askInt("Enter an Integer with Default value", 14.);
  message("Value read = %d\n",ianswer);
  
  ranswer = askDouble("Enter a Double with no Default value");
  message("Value read = %lf\n",ranswer);

  ranswer = askDouble("Enter a Double with Default value", 14.);
  message("Value read = %lf\n",ranswer);


  return 0;
}
