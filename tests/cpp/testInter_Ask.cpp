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
#include "Basic/File.hpp"
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

  // Testing string conversion Old to New style functions

  char buf_char[] = "abcde";
  message("Buffer (char) = %s\n",buf_char);
  message("Buffer size (adding the terminaison value) = %d\n",gslArraySize(buf_char));
  String buf_string("abcde");
  message("Buffer (string) = %s\n",buf_string.c_str());
  message("Buffer size = %d\n",static_cast<int>(buf_string.size()));

  // Testing strcpy

  message("Testing gslStrcpy:\n");
  message("- Expected result: abcde\n");

  char dst_char[10];
  gslStrcpy(dst_char,gslArraySize(dst_char),buf_char);
  message("- in char = %s\n",dst_char);

  String dst_string("fghijk");
  gslStrcpy(dst_string, buf_string);
  message("- in string = %s\n",dst_string.c_str());

  // Testing strcat

  message("Testing gslStrcat:\n");
  message("- Expected result: abcdefgh\n");

  char newbuf_char[] = "fgh";
  gslStrcat(dst_char,gslArraySize(dst_char),newbuf_char);
  message("- in char = %s\n",dst_char);

  String newbuf_string("fgh");
  gslStrcat(dst_string, newbuf_string);
  message("- in string = %s\n",dst_string.c_str());

  // Testing sprintf

  message("Testing gslSPrintf:\n");
  message("- Expected result: A+fgh\n");

  gslSPrintf(dst_char,gslArraySize(dst_char),"A+%s",newbuf_char);
  message("- in char = %s\n",dst_char);

  gslSPrintf(dst_string, "A+%s", newbuf_char);
  message("- in string = %s\n",dst_string.c_str());

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
