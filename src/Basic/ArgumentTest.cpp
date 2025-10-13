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
#include "Basic/ArgumentTest.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorNumT.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSquare.hpp"

ENUM_DEFINE(ENUM_TESTS)

namespace gstlrn
{

void _endOfLine()
{
  message("\n");
}
void _test()
{
  message("NA ");
}
void _introduction(const String& title, bool end_of_line = false)
{
  message("Testing for %s : ", title.c_str());
  if (end_of_line) _endOfLine();
}

void _printEmpty()
{
  message("Found an empty argument. This is correct\n");
}

void _printInt(Id value)
{
  if (IFFFF(value))
    _test();
  else
    message("%ld ", value);
}

void _printDouble(double value)
{
  if (FFFF(value))
    _test();
  else
    message("%lf ", value);
}

void _printString(const String& value)
{
  message("%s ", value.c_str());
}

void _printVectorInt(const VectorInt& values)
{
  for (Id i = 0; i < static_cast<Id>(values.size()); i++)
    _printInt(values[i]);
}

void _printVectorDouble(const VectorDouble& values)
{
  for (Id i = 0; i < static_cast<Id>(values.size()); i++)
    _printDouble(values[i]);
}

void _printVectorString(const VectorString& values)
{
  for (Id i = 0; i < static_cast<Id>(values.size()); i++)
    _printString(values[i]);
}

void _printVectorVectorInt(const VectorVectorInt& values)
{
  for (Id i = 0; i < static_cast<Id>(values.size()); i++)
  {
    for (Id j = 0; j < static_cast<Id>(values[i].size()); j++)
    {
      message("[%d][%d] : ", j + 1, i + 1);
      _printInt(values[i][j]);
      _endOfLine();
    }
  }
}

void _printVectorVectorDouble(const VectorVectorDouble& values)
{
  for (Id i = 0; i < static_cast<Id>(values.size()); i++)
  {
    for (Id j = 0; j < static_cast<Id>(values[i].size()); j++)
    {
      message("[%d][%d] : ", j + 1, i + 1);
      _printDouble(values[i][j]);
      _endOfLine();
    }
  }
}

/**
 * Function to test Integer argument
 * @param value Integer input argument
 */
void argumentTestInt(Id value)
{
  _introduction("TestInt");
  _printInt(value);
  _endOfLine();
}

/**
 * Function to test Double argument
 * @param value Double input argument
 */
void argumentTestDouble(double value)
{
  _introduction("TestDouble");
  _printDouble(value);
  _endOfLine();
}

void argumentTestVectorInt(const VectorInt& values)
{
  _introduction("TestVectorInt");
  _printVectorInt(values);
  _endOfLine();
}

void argumentTestVectorDouble(const VectorDouble& values)
{
  _introduction("TestVectorDouble");
  _printVectorDouble(values);
  _endOfLine();
}

void argumentTestString(const String& value)
{
  _introduction("TestString");
  _printString(value);
  _endOfLine();
}

void argumentTestVectorVectorInt(const VectorVectorInt& values)
{
  _introduction("TestVectorVectorInt", true);
  _printVectorVectorInt(values);
  _endOfLine();
}

void argumentTestVectorVectorDouble(const VectorVectorDouble& values)
{
  _introduction("TestVectorVectorDouble", true);
  _printVectorVectorDouble(values);
  _endOfLine();
}

void argumentTestVectorString(const VectorString& values)
{
  _introduction("TestVectorString", true);
  _printVectorString(values);
  _endOfLine();
}

void argumentTestStringOverload(const String& value)
{
  _introduction("TestString (Overload)");
  _printString(value);
  _endOfLine();
}

void argumentTestIntOverload(Id value)
{
  _introduction("TestInt (Overload)");
  _printInt(value);
  _endOfLine();
}

void argumentTestIntOverload(const VectorInt& values)
{
  _introduction("TestVectorInt (Overload)");
  _printVectorInt(values);
  _endOfLine();
}

void argumentTestDoubleOverload(double value)
{
  _introduction("TestDouble (Overload)");
  _printDouble(value);
  _endOfLine();
}

void argumentTestDoubleOverload(const VectorDouble& values)
{
  _introduction("TestVectorDouble (Overload)");
  _printVectorDouble(values);
  _endOfLine();
}

void argumentTestStringOverload(const VectorString& values)
{
  _introduction("TestVectorString (Overload)");
  _printVectorString(values);
  _endOfLine();
}

void argumentTestEnum(const ETests& value)
{
  message("Case : Value = %d - Descr = %s\n", value.getValue(), value.getDescr().data());
}

Id argumentReturnInt(Id value)
{
  _introduction("ReturnInt");
  _printInt(value);
  _endOfLine();
  return value;
}

double argumentReturnDouble(double value)
{
  _introduction("ReturnDouble");
  _printDouble(value);
  _endOfLine();
  return value;
}

VectorInt argumentReturnVectorInt(const VectorInt& values)
{
  _introduction("ReturnVectorInt");
  _printVectorInt(values);
  _endOfLine();
  return values;
}

VectorDouble argumentReturnVectorDouble(const VectorDouble& values)
{
  _introduction("ReturnVectorDouble");
  _printVectorDouble(values);
  _endOfLine();
  return values;
}

VectorVectorInt argumentReturnVectorVectorInt(const VectorVectorInt& values)
{
  _introduction("ReturnVectorVectorInt", true);
  _printVectorVectorInt(values);
  _endOfLine();
  return values;
}

VectorVectorDouble argumentReturnVectorVectorDouble(const VectorVectorDouble& values)
{
  _introduction("ReturnVectorVectorDouble", true);
  _printVectorVectorDouble(values);
  _endOfLine();
  return values;
}

void argumentDefTestInt(Id argInt)
{
  _introduction("DefTestInt");
  _printInt(argInt);
  _endOfLine();
}

void argumentDefTestDbl(double argDbl)
{
  _introduction("DefTestDouble");
  _printDouble(argDbl);
  _endOfLine();
}

void argumentDefTestStr(const String& argstr)
{
  _introduction("DefTestString");
  _printString(argstr);
  _endOfLine();
}

void argumentDefTestVInt(const VectorInt& argVInt)
{
  _introduction("DefTestVectorInt");
  if (argVInt.empty()) _printEmpty();
}

void argumentDefTestVDbl(const VectorDouble& argVDbl)
{
  _introduction("DefTestVectorDouble");
  if (argVDbl.empty()) _printEmpty();
}

void argumentDefTestVString(const VectorString& argVString)
{
  _introduction("DefTestVectorString");
  if (argVString.empty()) _printEmpty();
}

void argumentDefTestVVDbl(VectorVectorDouble argVVDbl)
{
  _introduction("DefTestVectorVectorDouble");
  if (argVVDbl.empty() || argVVDbl[0].empty()) _printEmpty();
}

void argumentDefTestVVInt(VectorVectorInt argVVInt)
{
  _introduction("DefTestVectorVectorInt");
  if (argVVInt.empty() || argVVInt[0].empty()) _printEmpty();
}

void argumentTestMatrixDense(const MatrixDense& mat)
{
  if (!mat.empty()) mat.display();
}
void argumentTestMatrixSquare(const MatrixSquare& mat)
{
  if (!mat.empty()) mat.display();
}
void argumentTestMatrixSymmetric(const MatrixSymmetric& mat)
{
  if (!mat.empty()) mat.display();
}
MatrixDense argumentReturnMatrix(Id nrows,
                                 Id ncols,
                                 Id seed)
{
  MatrixDense mat(nrows, ncols);
  mat.fillRandom(seed);
  return mat;
}

void argumentTestMatrixSparse(const MatrixSparse& mat)
{
  if (!mat.empty()) mat.display();
}

MatrixSparse argumentReturnMatrixSparse(Id nrows,
                                        Id ncols,
                                        double zeroPercent,
                                        Id seed)
{
  MatrixSparse mat(nrows, ncols);
  mat.fillRandom(seed, zeroPercent);
  mat.display();
  return mat;
}
} // namespace gstlrn