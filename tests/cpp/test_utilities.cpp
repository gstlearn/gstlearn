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
#include "geoslib_old_f.h"

#include "Basic/Utilities.hpp"

/****************************************************************************/
/*!
** Main Program for testing the utilities
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Utilities-");

  // =======================
  // Checking some functions
  // =======================

  double a = 1234.5678;
  double b = 0.000012345678;
  message("Value of a = %20.13lf\n",a);
  message("Value of b = %20.13lf\n",b);

  // Use of round operations

  mestitle(1, "RoundDecimals (all numbers are printed with 10 decimals)");
  for (int i = -3; i <= 4; i++)
    message("roundDecimals( a,%2d) = %20.10lf\n", i, truncateDecimals(a,i));
  for (int i = -3; i <= 4; i++)
    message("roundDecimals(-a,%2d) = %20.10lf\n", i, truncateDecimals(-a,i));
  for (int i = 5; i <= 10; i++)
    message("roundDecimals( b,%2d) = %20.10lf\n", i, truncateDecimals(b,i));
  for (int i = 5; i <= 10; i++)
    message("roundDecimals(-b,%2d) = %20.10lf\n", i, truncateDecimals(-b,i));

  mestitle(1, "RoundDigits (all numbers are printed with 13 decimals)");
  for (int i = 1; i <= 8; i++)
    message("roundDigits( a,%d) = %20.13lf\n", i, truncateDigits(a,i));
  for (int i = 1; i <= 8; i++)
    message("roundDigits(-a,%d) = %20.13lf\n", i, truncateDigits(-a,i));
  for (int i = 1; i <= 8; i++)
    message("roundDigits( b,%d) = %20.13lf\n", i, truncateDigits(b,i));
  for (int i = 1; i <= 8; i++)
    message("roundDigits(-b,%d) = %20.13lf\n", i, truncateDigits(-b,i));

  return(0);
}
