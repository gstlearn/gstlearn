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
void _test()
{
  message("NA ");
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

void argumentTestInt(int value)
{
  _introduction("Integer");
  _printInt(value);
  _endOfLine();
}

void argumentTestDouble(double value)
{
  _introduction("Double");
  _printDouble(value);
  _endOfLine();
}

void argumentTestVectorInt(const VectorInt& values)
{
  _introduction("Vector Integer");
  for (int i = 0; i < (int) values.size(); i++)
    _printInt(values[i]);
  _endOfLine();
}

void argumentTestVectorDouble(const VectorDouble& values)
{
  _introduction("Vector Double");
  for (int i = 0; i < (int) values.size(); i++)
    _printDouble(values[i]);
  _endOfLine();
}

void argumentTestString(const String& value)
{
  _introduction("String");
  _printString(value);
  _endOfLine();
}

void argumentTestVectorString(const VectorString& values)
{
  _introduction("Vector String");
  for (int i = 0; i < (int) values.size(); i++)
    _printString(values[i]);
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
  _introduction("Vector Integer (Overload)");
  for (int i = 0; i < (int) values.size(); i++)
    _printInt(values[i]);
  _endOfLine();
}

void argumentTestStringOverload(const VectorString& values)
{
  _introduction("Vector String (Overload)");
  for (int i = 0; i < (int) values.size(); i++)
    _printString(values[i]);
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

