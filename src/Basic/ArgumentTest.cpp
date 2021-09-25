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
  std::cout << "Testing for " << title << " : ";
}
void _printInt(int value)
{
  if (IFFFF(value))
    std::cout << "NA" << " ";
  else
    std::cout << value << " ";
}
void _printDouble(double value)
{
  if (FFFF(value))
    std::cout << "NA" << " ";
  else
    std::cout << value << " ";
}
void _printString(const String& value)
{
  std::cout << value << " ";
}
void _endOfLine()
{
  std::cout << std::endl;
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

void argumentTestSurcharge(const String& value)
{
  _introduction("String (surcharge)");
  _printString(value);
  _endOfLine();
}

void argumentTestSurcharge(const VectorString& values)
{
  _introduction("Vector String (surcharge)");
  for (int i = 0; i < (int) values.size(); i++)
    _printString(values[i]);
  _endOfLine();
}

void argumentTestEnum(ETests value)
{
  std::cout << value << std::endl; // Automatic conversion to int?
}

