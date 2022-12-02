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
#include "Style.hpp"

// For message function
#include "Basic/AStringable.hpp"

/**
 * Default constructor
 *
 * Always use the initialization list to set the default value of each class member.
 * Put attributes in the same order than their declaration in the C++ header.
 */
Style::Style()
    : _argInt(0),
      _argDouble(0.),
      _argVectorInt(),
      _argVectorDouble()
{
}

/**
 * Copy constructor
 *
 * Always use the initialization list to set the value of each class member.
 * Put attributes in the same order than their declaration in the C++ header.
 */
Style::Style(const Style& r)
    : _argInt(r._argInt),
      _argDouble(r._argDouble),
      _argVectorInt(r._argVectorInt),
      _argVectorDouble(r._argVectorDouble)
{
}

/**
 * Assignment operator
 *
 * Always protect from self copy if the if statement
 */
Style& Style::operator=(const Style &r)
{
  if (this != &r)
  {
    _argInt = r._argInt;
    _argDouble = r._argDouble;
    _argVectorInt = r._argVectorInt;
    _argVectorDouble = r._argVectorDouble;
  }
  return *this;
}

/**
 * Destructor
 *
 * Always free pointers and clear the lists of the class
 */
Style::~Style()
{
}

/**
 * A method with standard argument documentation
 *
 * @param myArg Here should be placed the description of this argument
 *
 * @return Description of the returned value
 */
int Style::DocumentedStandard(int myArg) const
{
  message("Documented Function\n");
  message("Value of MyArg = %d\n",myArg);
  return 0;
}

/**
 * Documentation with Latex formula
 *
 * The distance between \f$ p1=(x_1,y_1) \f$ and \f$ p2=(x_2,y_2) \f$ is \f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}\f$
 *
 * @see SpaceRN::getDistance
 *
 * @param myArg Here should be placed the description of this argument
 *
 * @return The value of the argument + 1
 */
int Style::DocumentedWithFormula(int myArg) const
{
  message("Input Argument = %d\n", myArg);
  return _increment(myArg);
}

/**
 * Description of a static function
 *
 * @param myArgInt Integer argument
 * @param myArgDoubleDef Double argument
 */
void Style::myFunction(int myArgInt, double myArgDoubleDef)
{
  message("Input Argument (int) = %d\n", myArgInt);
  message("Input Argument (double with default value) = %lf\n", myArgDoubleDef);
}

/**
 * A function where the argument is not used (could be the case in abstract methods)
 *
 * @param a Input argument not used (but documented)
 *
 * @return Error returned code
 */
int Style::UnusedArgument(int a)
{
  // Use the SYMBO_UNSED macro to prevent compiler warning
  // Do not comment the argument (NOT /*a*/)
  SYMBOL_UNUSED(a);
  return 0;
}

/**
 *  Example of a private method
 *
 *  @param arg Input integer argument
 *
 *  @return The input argument incremented
 */
int Style::_increment(int arg) const
{
  return arg++;
}
