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

#include "Basic/Vector.hpp"

bool isInteger(double value, double eps = EPSILON10);
int  getClosestInteger(double value);
bool isMultiple(int nbig, int nsmall);
bool isOdd(int number);
bool isEven(int number);
int  FFFF(double value);
int  IFFFF(int value);
double getTEST();
int    getITEST();
double getMin(double val1, double val2);
double getMax(double val1, double val2);
