/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
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
