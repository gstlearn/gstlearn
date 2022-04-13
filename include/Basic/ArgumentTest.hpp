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
#include "Basic/Vector.hpp"
#include "Enum/AEnum.hpp"

#define ENUM_TESTS ETests, CASE0,\
                   CASE0 , 0, "Enum test case 0",\
                   CASE1 , 1, "Enum test case 1",\
                   CASE2 , 2, "Enum test case 2"

ENUM_DECLARE(ENUM_TESTS)

class GSTLEARN_EXPORT argClass
{
public:
  int    _ival;
  double _rval;
  String _sval;

  argClass(int ival = -1, double rval = -1.1, String sval = "NA")
      : _ival(ival),
        _rval(rval),
        _sval(sval)
  {
  }

  int    getIval() const { return _ival; }
  void   setIval(int ival) { _ival = ival; }
  double getRval() const { return _rval; }
  void   setRval(double rval) { _rval = rval; }
  const String& getSval() const { return _sval; }
  void   setSval(const String& sval) { _sval = sval; }
  void   display() const {
    message("Integer = %d - Real = %lf - String = %s\n", _ival,_rval,_sval.c_str());
  }
};

/**
 * This file is meant to provide a set of functions for testing arguments
 * for Python and R interfaces
 */

GSTLEARN_EXPORT void argumentTestInt(int value);
GSTLEARN_EXPORT void argumentTestDouble(double value);
GSTLEARN_EXPORT void argumentTestVectorInt(const VectorInt& values);
GSTLEARN_EXPORT void argumentTestVectorDouble(const VectorDouble& values);
GSTLEARN_EXPORT void argumentTestString(const String& value);
GSTLEARN_EXPORT void argumentTestVectorString(const VectorString& values);

GSTLEARN_EXPORT void argumentTestIntOverload(int value);
GSTLEARN_EXPORT void argumentTestIntOverload(const VectorInt& values);
GSTLEARN_EXPORT void argumentTestStringOverload(const String& value);
GSTLEARN_EXPORT void argumentTestStringOverload(const VectorString& values);

GSTLEARN_EXPORT void argumentTestEnum(ETests value);

GSTLEARN_EXPORT int argumentReturnInt(int value);
GSTLEARN_EXPORT double argumentReturnDouble(double value);
