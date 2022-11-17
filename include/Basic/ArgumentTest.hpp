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
#include "Enum/AEnum.hpp"
#include "Basic/Utilities.hpp"

#define ENUM_TESTS ETests, CASE0,\
                   CASE0 , 0, "Enum test case 0",\
                   CASE1 , 1, "Enum test case 1",\
                   CASE2 , 2, "Enum test case 2"

ENUM_DECLARE(ENUM_TESTS)

class GSTLEARN_EXPORT argClass
{
  // The members are public (just for testing)
public:
  int    ival;
  double rval;
  String sval;

  argClass(int iival = -1, double rrval = -1.1, String ssval = STRING_NA)
      : ival(iival),
        rval(rrval),
        sval(ssval)
  {
  }

  int    getIval() const { return ival; }
  void   setIval(int iival) { ival = iival; }
  double getRval() const { return rval; }
  void   setRval(double rrval) { rval = rrval; }
  const String& getSval() const { return sval; }
  void   setSval(const String& ssval) { sval = ssval; }
  void   display() const {
    message("Integer = %d - Real = %lf - String = %s\n", ival, rval, sval.c_str());
  }
};

GSTLEARN_EXPORT void argumentTestInt(int value);
GSTLEARN_EXPORT void argumentTestDouble(double value);
GSTLEARN_EXPORT void argumentTestVectorInt(const VectorInt& values);
GSTLEARN_EXPORT void argumentTestVectorDouble(const VectorDouble& values);
GSTLEARN_EXPORT void argumentTestVectorVectorInt(const VectorVectorInt& values);
GSTLEARN_EXPORT void argumentTestVectorVectorDouble(const VectorVectorDouble& values);
GSTLEARN_EXPORT void argumentTestString(const String& value);
GSTLEARN_EXPORT void argumentTestVectorString(const VectorString& values);

GSTLEARN_EXPORT void argumentTestIntOverload(int value);
GSTLEARN_EXPORT void argumentTestIntOverload(const VectorInt& values);
GSTLEARN_EXPORT void argumentTestDoubleOverload(double value);
GSTLEARN_EXPORT void argumentTestDoubleOverload(const VectorDouble& values);
GSTLEARN_EXPORT void argumentTestStringOverload(const String& value);
GSTLEARN_EXPORT void argumentTestStringOverload(const VectorString& values);

GSTLEARN_EXPORT void argumentTestEnum(ETests value);

GSTLEARN_EXPORT int argumentReturnInt(int value);
GSTLEARN_EXPORT double argumentReturnDouble(double value);
GSTLEARN_EXPORT VectorInt argumentReturnVectorInt(const VectorInt& values);
GSTLEARN_EXPORT VectorDouble argumentReturnVectorDouble(const VectorDouble& values);

GSTLEARN_EXPORT void argumentDefTestInt(int argInt = 2);
GSTLEARN_EXPORT void argumentDefTestDbl(double argDbl = 2.);
GSTLEARN_EXPORT void argumentDefTestStr(String argstr = "Default String");
GSTLEARN_EXPORT void argumentDefTestVInt(VectorInt argVInt = VectorInt());
GSTLEARN_EXPORT void argumentDefTestVDbl(VectorDouble argdVDbl = VectorDouble());
GSTLEARN_EXPORT void argumentDefTestVString(VectorString argVString = VectorString());
GSTLEARN_EXPORT void argumentDefTestVVInt(VectorVectorInt argVVInt = VectorVectorInt());
GSTLEARN_EXPORT void argumentDefTestVVDbl(VectorVectorDouble argVVDbl = VectorVectorDouble());
