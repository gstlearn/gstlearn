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

/**
 * This test is meant to check the interactive questioning
 * @return
 */
int main()
{
  // This (interactive) test cannot be performed in batch
  bool flagPerform = false;
  if (! flagPerform)
  {
    message("This test is not performed on purpose\n");
    message("as it does not make sense in Batch\n");
    return 0;
  }

  int ianswer;
  double ranswer;

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
