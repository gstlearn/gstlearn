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
#include <cmath>
#include "geoslib_d.h"
#include "Basic/Utilities.hpp"

GSTLEARN_EXPORT bool isInteger(double value, double eps)
{
  int iclose = getClosestInteger(value);
  if (ABS((double) iclose - value) > eps) return false;
  return true;
}

GSTLEARN_EXPORT int getClosestInteger(double value)
{
  int iclose = (int) round(value);
  return iclose;
}

GSTLEARN_EXPORT bool isMultiple(int nbig, int nsmall)
{
  double ratio;

  ratio = (double) nbig / (double) nsmall;
  return (isInteger(ratio));
}

GSTLEARN_EXPORT bool isOdd(int number)
{
  int middle;

  middle = number / 2;
  if (number != 2 * middle)
    return true;
  else
    return false;
}

GSTLEARN_EXPORT bool isEven(int number)
{
  int middle;

  middle = number / 2;
  if (number != 2 * middle)
    return false;
  else
    return true;
}

GSTLEARN_EXPORT double getMin(double val1, double val2)
{
  if (FFFF(val1)) return (val2);
  if (FFFF(val2)) return (val1);
  return (MIN(val1, val2));
}

GSTLEARN_EXPORT double getMax(double val1, double val2)
{
  if (FFFF(val1)) return (val2);
  if (FFFF(val2)) return (val1);
  return (MAX(val1, val2));
}

GSTLEARN_EXPORT double getTEST()
{
  return TEST;
}

GSTLEARN_EXPORT int    getITEST()
{
  return ITEST;
}

/****************************************************************************/
/*!
 **  Checks if a double value is TEST
 **
 ** \return  1 if a TEST value is encountered; 0 otherwise
 **
 ** \param[in]  value Value to be tested
 **
 *****************************************************************************/
GSTLEARN_EXPORT int FFFF(double value)
{
  int rep;

  rep = 0;
  if (std::isnan(value)) rep = 1;
  if (value > TEST_COMP) rep = 1;

  return (rep);
}

/****************************************************************************/
/*!
 **  Checks if an integer value is TEST
 **
 ** \return  1 if a ITEST value is encountered; 0 otherwise
 **
 ** \param[in]  value Value to be tested
 **
 *****************************************************************************/
GSTLEARN_EXPORT int IFFFF(int value)
{
  int rep;

  rep = 0;
  if (value == ITEST) rep = 1;

  return (rep);
}
