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
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/AStringable.hpp"
#include <iostream>

/**
 * This test is meant to check the string manipulations
 */
int main()
{
  // Testing string conversion Old to New style functions

  char buf_char[] = "abcde";
  message("Buffer (char) = %s\n",buf_char);

  // Testing strcpy

  message("Testing gslStrcpy:\n");
  message("- Expected result: abcde\n");

  char dst_char[10];
  gslStrcpy(dst_char,buf_char);
  message("- Result = %s\n",dst_char);

  // Testing strcat

  message("Testing gslStrcat:\n");
  message("- Expected result: abcdefgh\n");

  char newbuf_char[] = "fgh";
  gslStrcat(dst_char,newbuf_char);
  message("- Result = %s\n",dst_char);

  // Testing sprintf

  message("Testing gslSPrintf:\n");
  message("- Expected result: A+fgh\n");
  gslSPrintf(dst_char,"A+%s",newbuf_char);
  message("- Result(char) = %s\n",dst_char);

  return 0;
}
