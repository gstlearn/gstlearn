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
#pragma once

#include "gstlearn_export.hpp"
#include "Enum/AEnum.hpp"
#include "Basic/Utilities.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixSparse.hpp"

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

  argClass(int iival = -1, double rrval = -1.1, const String& ssval = STRING_NA)
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

GSTLEARN_EXPORT void argumentTestEnum(const ETests& value);

GSTLEARN_EXPORT int argumentReturnInt(int value);
GSTLEARN_EXPORT double argumentReturnDouble(double value);
GSTLEARN_EXPORT VectorInt argumentReturnVectorInt(const VectorInt& values);
GSTLEARN_EXPORT VectorDouble argumentReturnVectorDouble(const VectorDouble& values);
GSTLEARN_EXPORT VectorVectorInt argumentReturnVectorVectorInt(const VectorVectorInt& values);
GSTLEARN_EXPORT VectorVectorDouble argumentReturnVectorVectorDouble(const VectorVectorDouble& values);

GSTLEARN_EXPORT void argumentDefTestInt(int argInt = 2);
GSTLEARN_EXPORT void argumentDefTestDbl(double argDbl = 2.);
GSTLEARN_EXPORT void argumentDefTestStr(const String& argstr = "Default String");
GSTLEARN_EXPORT void argumentDefTestVInt(const VectorInt& argVInt = VectorInt());
GSTLEARN_EXPORT void argumentDefTestVDbl(const VectorDouble& argVDbl = VectorDouble());
GSTLEARN_EXPORT void argumentDefTestVString(const VectorString& argVString = VectorString());
GSTLEARN_EXPORT void argumentDefTestVVInt(VectorVectorInt argVVInt = VectorVectorInt());
GSTLEARN_EXPORT void argumentDefTestVVDbl(VectorVectorDouble argVVDbl = VectorVectorDouble());

GSTLEARN_EXPORT void argumentTestMatrixDense(const MatrixDense& mat = MatrixDense());
GSTLEARN_EXPORT void argumentTestMatrixSquareGeneral(const MatrixSquareGeneral& mat = MatrixSquareGeneral());
GSTLEARN_EXPORT void argumentTestMatrixSymmetric(const MatrixSymmetric& mat = MatrixSymmetric());

GSTLEARN_EXPORT MatrixDense argumentReturnMatrix(int nrows = 2,
                                                 int ncols = 3,
                                                 int seed  = 1312);

GSTLEARN_EXPORT void argumentTestMatrixSparse(const MatrixSparse& mat = MatrixSparse());
GSTLEARN_EXPORT MatrixSparse argumentReturnMatrixSparse(int nrows = 2,
                                                        int ncols = 3,
                                                        double zeroPercent = 0.1,
                                                        int seed = 1356);
