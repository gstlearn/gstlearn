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

typedef enum
{
  TEST_CASE0 = 0,
  TEST_CASE1 = 1,
  TEST_CASE2 = 2,
} ENUM_TESTS;

/**
 * This file is meant to provide a set of functions for testing arguments
 * for Python and R interfaces
 */

void argumentTestInt(int value);
void argumentTestDouble(double value);
void argumentTestVectorInt(const VectorInt& values);
void argumentTestVectorDouble(const VectorDouble& values);
void argumentTestString(const String& value);
void argumentTestVectorString(const VectorString& values);

void argumentTestSurcharge(const String& value);
void argumentTestSurcharge(const VectorString& values);

void argumentTestEnum(ENUM_TESTS value);
