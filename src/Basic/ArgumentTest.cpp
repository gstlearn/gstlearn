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
#include "Matrix/AMatrix.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/ArgumentTest.hpp"

ENUM_DEFINE(ENUM_TESTS)

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
  message("Testing for %s : ",title.c_str());
  if (end_of_line) _endOfLine();
}

void _printEmpty()
{
  message("Found an empty argument. This is correct\n");
}

void _printInt(int value)
{
  if (IFFFF(value))
   _test();
  else
    message("%d ",value);
}

void _printDouble(double value)
{
  if (FFFF(value))
    _test();
  else
    message("%lf ",value);
}

void _printString(const String& value)
{
  message("%s ",value.c_str());
}

void _printVectorInt(const VectorInt& values)
{
  for (int i = 0; i < (int) values.size(); i++)
    _printInt(values[i]);
}

void _printVectorDouble(const VectorDouble& values)
{
  for (int i = 0; i < (int) values.size(); i++)
    _printDouble(values[i]);
}

void _printVectorString(const VectorString& values)
{
  for (int i = 0; i < (int) values.size(); i++)
    _printString(values[i]);
}

void _printVectorVectorInt(const VectorVectorInt& values)
{
  for (int i = 0; i < (int) values.size(); i++)
  {
    for (int j = 0; j < (int) values[i].size(); j++)
    {
      message("[%d][%d] : ",j+1,i+1);
      _printInt(values[i][j]);
      _endOfLine();
    }
  }
}

void _printVectorVectorDouble(const VectorVectorDouble& values)
{
  for (int i = 0; i < (int) values.size(); i++)
  {
    for (int j = 0; j < (int) values[i].size(); j++)
    {
      message("[%d][%d] : ",j+1,i+1);
      _printDouble(values[i][j]);
      _endOfLine();
    }
  }
}

/**
 * Function to test Integer argument
 * @param value Integer input argument
 */
void argumentTestInt(int value)
{
  _introduction("Integer");
  _printInt(value);
  _endOfLine();
}

/**
 * Function to test Double argument
 * @param value Double input argument
 */
void argumentTestDouble(double value)
{
  _introduction("Double");
  _printDouble(value);
  _endOfLine();
}

void argumentTestVectorInt(const VectorInt& values)
{
  _introduction("VectorInt");
  _printVectorInt(values);
  _endOfLine();
}

void argumentTestVectorDouble(const VectorDouble& values)
{
  _introduction("VectorDouble");
  _printVectorDouble(values);
  _endOfLine();
}

void argumentTestString(const String& value)
{
  _introduction("String");
  _printString(value);
  _endOfLine();
}

void argumentTestVectorVectorInt(const VectorVectorInt& values)
{
  _introduction("VectorVectorInt",true);
  _printVectorVectorInt(values);
  _endOfLine();
}

void argumentTestVectorVectorDouble(const VectorVectorDouble& values)
{
  _introduction("VectorVectorDouble",true);
  _printVectorVectorDouble(values);
  _endOfLine();
}

void argumentTestVectorString(const VectorString& values)
{
  _introduction("VectorString");
  _printVectorString(values);
  _endOfLine();
}

void argumentTestStringOverload(const String& value)
{
  _introduction("String (Overload)");
  _printString(value);
  _endOfLine();
}

void argumentTestIntOverload(int value)
{
  _introduction("Int (Overload)");
  _printInt(value);
  _endOfLine();
}

void argumentTestIntOverload(const VectorInt& values)
{
  _introduction("VectorInt (Overload)");
  _printVectorInt(values);
  _endOfLine();
}

void argumentTestDoubleOverload(double value)
{
  _introduction("Double (Overload)");
  _printDouble(value);
  _endOfLine();
}

void argumentTestDoubleOverload(const VectorDouble& values)
{
  _introduction("VectorDouble (Overload)");
  _printVectorDouble(values);
  _endOfLine();
}

void argumentTestStringOverload(const VectorString& values)
{
  _introduction("VectorString (Overload)");
  _printVectorString(values);
  _endOfLine();
}

void argumentTestEnum(const ETests& value)
{
  message("Case : Value = %d - Descr = %s\n", value.getValue(),value.getDescr().c_str());
}

int argumentReturnInt(int value)
{
  _introduction("Integer");
  _printInt(value);
  _endOfLine();
  return value;
}

double argumentReturnDouble(double value)
{
  _introduction("Double");
  _printDouble(value);
  _endOfLine();
  return value;
}

VectorInt argumentReturnVectorInt(const VectorInt& values)
{
  _introduction("VectorInt");
  _printVectorInt(values);
  _endOfLine();
  return values;
}

GSTLEARN_EXPORT VectorDouble argumentReturnVectorDouble(const VectorDouble& values)
{
  _introduction("VectorDouble");
  _printVectorDouble(values);
  _endOfLine();
  return values;
}

GSTLEARN_EXPORT void argumentDefTestInt(int argInt)
{
  _introduction("Integer");
  _printInt(argInt);
  _endOfLine();
}

GSTLEARN_EXPORT void argumentDefTestDbl(double argDbl)
{
  _introduction("Double");
   _printDouble(argDbl);
   _endOfLine();
}

GSTLEARN_EXPORT void argumentDefTestStr(const String& argstr)
{
  _introduction("String");
   _printString(argstr);
   _endOfLine();
}

GSTLEARN_EXPORT void argumentDefTestVInt(const VectorInt& argVInt)
{
  _introduction("Vector Int");
  if (argVInt.empty()) _printEmpty();
}

GSTLEARN_EXPORT void argumentDefTestVDbl(const VectorDouble& argVDbl)
{
  _introduction("Vector Double");
  if (argVDbl.empty()) _printEmpty();
}

GSTLEARN_EXPORT void argumentDefTestVString(const VectorString& argVString)
{
  _introduction("Vector String");
  if (argVString.empty()) _printEmpty();
}

GSTLEARN_EXPORT void argumentDefTestVVDbl(VectorVectorDouble argVVDbl)
{
  _introduction("Vector Vector Double");
  if (argVVDbl.empty() || argVVDbl[0].empty()) _printEmpty();
}

GSTLEARN_EXPORT void argumentDefTestVVInt(VectorVectorInt argVVInt)
{
  _introduction("Vector Vector Int");
  if (argVVInt.empty() || argVVInt[0].empty()) _printEmpty();
}

GSTLEARN_EXPORT void argumentTestMatrixRectangular(const MatrixRectangular& mat)
{
  if (mat.empty())
    _printEmpty();
  else
    mat.display();
}
