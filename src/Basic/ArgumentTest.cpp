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
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Vector.hpp"
#include "Basic/ArgumentTest.hpp"

ENUM_DEFINE(ENUM_TESTS)

void _introduction(const String& title)
{
  message("Testing for %s : ",title.c_str());
}
void _endOfLine()
{
  message("\n");
}
void _nextOfLine()
{
  message(" - ");
}
void _test()
{
  message("NA ");
}

void _printInt(int value, bool next = false)
{
  if (IFFFF(value))
   _test();
  else
    message("%d",value);
  if (next) _nextOfLine();
}
void _printDouble(double value, bool next = false)
{
  if (FFFF(value))
    _test();
  else
    message("%lf",value);
  if (next) _nextOfLine();
}
void _printString(const String& value, bool next = false)
{
  message("%s",value.c_str());
  if (next) _nextOfLine();
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
  for (int i = 0; i < (int) values.size(); i++)
    _printInt(values[i],true);
  _endOfLine();
}

void argumentTestVectorDouble(const VectorDouble& values)
{
  _introduction("VectorDouble");
  for (int i = 0; i < (int) values.size(); i++)
    _printDouble(values[i],true);
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
  _introduction("VectorVectorInt");
  message("Dimension First Level = %d\n",(int) values.size());
  for (auto &e: values)
  {
    message("Dimension of Second Level = %d\n",(int) e.size());
    for (auto &f: e)
      _printInt(f, true);
    _endOfLine();
  }
}

void argumentTestVectorString(const VectorString& values)
{
  _introduction("VectorString");
  for (int i = 0; i < (int) values.size(); i++)
    _printString(values[i],true);
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
  for (int i = 0; i < (int) values.size(); i++)
    _printInt(values[i], true);
  _endOfLine();
}

void argumentTestStringOverload(const VectorString& values)
{
  _introduction("VectorString (Overload)");
  for (int i = 0; i < (int) values.size(); i++)
    _printString(values[i],true);
  _endOfLine();
}

void argumentTestEnum(ETests value)
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

