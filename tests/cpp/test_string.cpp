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
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/AStringable.hpp"

#include <iostream>

/**
 * This test is meant to check the string manipulations
 */
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

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
