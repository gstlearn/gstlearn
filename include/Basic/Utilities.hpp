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
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

GSTLEARN_EXPORT bool   isInteger(double value, double eps = EPSILON10);
GSTLEARN_EXPORT int    getClosestInteger(double value);
GSTLEARN_EXPORT bool   isMultiple(int nbig, int nsmall);
GSTLEARN_EXPORT bool   isOdd(int number);
GSTLEARN_EXPORT bool   isEven(int number);
GSTLEARN_EXPORT int    FFFF(double value);
GSTLEARN_EXPORT int    IFFFF(int value);
GSTLEARN_EXPORT double getTEST();
GSTLEARN_EXPORT int    getITEST();
GSTLEARN_EXPORT double getMin(double val1, double val2);
GSTLEARN_EXPORT double getMax(double val1, double val2);
